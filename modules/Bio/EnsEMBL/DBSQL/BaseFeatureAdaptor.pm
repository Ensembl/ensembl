#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor - Abstract Base class for 
                                          FeatureAdaptors

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for feature adaptors. This base class is simply a way
of eliminating code duplication through the implementation of methods 
common to all feature adaptors.

=head1 CONTACT

Contact EnsEMBL development list for info: <ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use vars qw(@ISA $SLICE_FEATURE_CACHE_SIZE);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Cache;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

$SLICE_FEATURE_CACHE_SIZE = 4;


=head2 new

  Arg [1]    : list of args @args
               Superclass constructor arguments
  Example    : none
  Description: Constructor which just initializes internal cache structures
  Returntype : Bio::EnsEMBL::BaseFeatureAdaptor
  Exceptions : none
  Caller     : implementing subclass constructors

=cut

sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(@args);

  #initialize caching data structures
  $self->{'_slice_feature_cache'} = {};

  tie(%{$self->{'_slice_feature_cache'}}, 
      'Bio::EnsEMBL::Utils::Cache',
      $SLICE_FEATURE_CACHE_SIZE);

  return $self;
}

=head2 generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) string $logic_name
               the logic_name of the analysis of the features to obtain
  Example    : $fts = $a->generic_fetch('contig_id in (1234, 1235)', 'Swall');
  Description: Performs a database fetch and returns feature objects in
               contig coordinates.
  Returntype : listref of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : none
  Caller     : BaseFeatureAdaptor, ProxyDnaAlignFeatureAdaptor::generic_fetch

=cut
  
sub generic_fetch {
  my ($self, $constraint, $logic_name, $mapper, $slice) = @_;
  
  my @tabs = $self->_tables;
  my $columns = join(', ', $self->_columns());
  
  if($logic_name) {
    #determine the analysis id via the logic_name
    my $analysis = 
      $self->db->get_AnalysisAdaptor()->fetch_by_logic_name($logic_name);
    unless(defined $analysis && $analysis->dbID() ) {
      $self->warn("No analysis for logic name $logic_name exists\n");
      return [];
    }
    
    my $analysis_id = $analysis->dbID();

    #get the synonym for the primary table
    my $syn = $tabs[0]->[1];

    if($constraint) {
      $constraint .= " AND ${syn}.analysis_id = $analysis_id";
    } else {
      $constraint = " ${syn}.analysis_id = $analysis_id";
    }
  } 

  #
  # Construct a left join statement if one was defined, and remove the
  # left-joined table from the table list
  #
  my ($tablename, $condition) = $self->_left_join;
  my $left_join = '';
  my @tables;
  if($tablename && $condition) {
    while(my $t = shift @tabs) {
      if($tablename eq $t->[0]) {
	my $syn = $t->[1]; 
	$left_join =  "LEFT JOIN $tablename $syn $condition";
	push @tables, @tabs;
	last;
      } else {
	push @tables, $t;
      }
    }
  } else {
    @tables = @tabs;
  }
      
  #construct a nice table string like 'table1 t1, table2 t2'
  my $tablenames = join(', ', map({ join(' ', @$_) } @tables));

  my $sql = "SELECT $columns FROM $tablenames $left_join";

  my $default_where = $self->_default_where_clause;
  my $final_clause = $self->_final_clause;

  #append a where clause if it was defined
  if($constraint) { 
    $sql .= " where $constraint ";
    if($default_where) {
      $sql .= " and $default_where ";
    }
  } elsif($default_where) {
    $sql .= " where $default_where ";
  }

  #append additional clauses which may have been defined
  $sql .= " $final_clause";

  my $sth = $self->prepare($sql);
  
  $sth->execute;  

  return $self->_objs_from_sth($sth, $mapper, $slice);
}


=head2 fetch_by_dbID

  Arg [1]    : int $id
               the unique database identifier for the feature to be obtained 
  Example    : $feat = $adaptor->fetch_by_dbID(1234);
  Description: Returns the feature created from the database defined by the
               the id $id. 
  Returntype : Bio::EnsEMBL::SeqFeature
  Exceptions : thrown if $id is not defined
  Caller     : general

=cut

sub fetch_by_dbID{
  my ($self,$id) = @_;
  
  unless(defined $id) {
    $self->throw("fetch_by_dbID must have an id");
  }

  my @tabs = $self->_tables;

  my ($name, $syn) = @{$tabs[0]};

  #construct a constraint like 't1.table1_id = 1'
  my $constraint = "${syn}.${name}_id = $id";

  #return first element of _generic_fetch list
  my ($feat) = @{$self->generic_fetch($constraint)}; 
  return $feat;
}


=head2 fetch_all_by_RawContig_constraint

  Arg [1]    : Bio::EnsEMBL::RawContig $contig
               The contig object from which features are to be obtained
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_Contig_constraint($ctg,'perc_ident>5.0');
  Description: Returns a listref of features created from the database which 
               are on the contig defined by $cid and fulfill the SQL constraint
               defined by $constraint. If logic name is defined, only features
               with an analysis of type $logic_name will be returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : thrown if $cid is not defined
  Caller     : general

=cut

sub fetch_all_by_RawContig_constraint {
  my ($self, $contig, $constraint, $logic_name) = @_;
  
  unless( defined $contig ) {
    $self->throw("fetch_by_Contig_constraint must have an contig");
  }

  unless( ref $contig && $contig->isa('Bio::EnsEMBL::RawContig')) {
    $self->throw("contig argument is not a Bio::EnsEMBL::RawContig object\n");
  }

  my $cid = $contig->dbID();

  #get the synonym of the primary_table
  my @tabs = $self->_tables;
  my $syn = $tabs[0]->[1];

  if($constraint) {
    $constraint .= " AND ${syn}.contig_id = $cid";
  } else {
    $constraint = "${syn}.contig_id = $cid";
  }

  return $self->generic_fetch($constraint, $logic_name);
}


=head2 fetch_all_by_RawContig

  Arg [1]    : Bio::EnsEMBL::RawContig $contig 
               the contig from which features should be obtained
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : @fts = $a->fetch_all_by_RawContig($contig, 'wall');
  Description: Returns a list of features created from the database which are 
               are on the contig defined by $cid If logic name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::*Feature in contig coordinates
  Exceptions : none
  Caller     : general

=cut
   
sub fetch_all_by_RawContig {
  my ( $self, $contig, $logic_name ) = @_;

  return $self->fetch_all_by_RawContig_constraint($contig, '',$logic_name);
}


=head2 fetch_all_by_RawContig_and_score
  Arg [1]    : Bio::EnsEMBL::RawContig $contig 
               the contig from which features should be obtained
  Arg [2]    : (optional) float $score
               the lower bound of the score of the features to obtain
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : @fts = $a->fetch_by_RawContig_and_score(1, 50.0, 'Swall');
  Description: Returns a list of features created from the database which are 
               are on the contig defined by $cid and which have score greater  
               than score.  If logic name is defined, only features with an 
               analysis of type $logic_name will be returned. 
  Returntype : listref of Bio::EnsEMBL::*Feature in contig coordinates
  Exceptions : thrown if $score is not defined
  Caller     : general

=cut

sub fetch_all_by_RawContig_and_score{
  my($self, $contig, $score, $logic_name) = @_;

  my $constraint;

  if(defined $score){
    #get the synonym of the primary_table
    my @tabs = $self->_tables;
    my $syn = $tabs[0]->[1];
    $constraint = "${syn}.score > $score";
  }
    
  return $self->fetch_all_by_RawContig_constraint($contig, $constraint, 
					       $logic_name);
}


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fts = $a->fetch_all_by_Slice($slice, 'Swall');
  Description: Returns a listref of features created from the database 
               which are on the Slice defined by $slice. If $logic_name is 
               defined only features with an analysis of type $logic_name 
               will be returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $logic_name) = @_;
  
  #fetch by constraint with empty constraint
  return $self->fetch_all_by_Slice_constraint($slice, '', $logic_name);
}


=head2 fetch_all_by_Slice_and_score

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) float $score
               lower bound of the the score of the features retrieved
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fts = $a->fetch_all_by_Slice($slice, 'Swall');
  Description: Returns a list of features created from the database which are 
               are on the Slice defined by $slice and which have a score 
               greated than $score. If $logic_name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice

=cut

sub fetch_all_by_Slice_and_score {
  my ($self, $slice, $score, $logic_name) = @_;
  my $constraint;

  if(defined $score) {
    #get the synonym of the primary_table
    my @tabs = $self->_tables;
    my $syn = $tabs[0]->[1];
    $constraint = "${syn}.score > $score";
  }

  return $self->fetch_all_by_Slice_constraint($slice, $constraint, 
					      $logic_name);
}  


=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'perc_ident > 5');
  Description: Returns a listref of features created from the database which 
               are on the Slice defined by $slice and fulfill the SQL 
               constraint defined by $constraint. If logic name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice

=cut

sub fetch_all_by_Slice_constraint {
  my($self, $slice, $constraint, $logic_name) = @_;

  unless(defined $slice && ref $slice && $slice->isa("Bio::EnsEMBL::Slice")) {
    $self->throw("Slice arg must be a Bio::EnsEMBL::Slice not a [$slice]\n");
  }

  $logic_name = '' unless $logic_name;
  $constraint = '' unless $constraint;

  #check the cache and return if we have already done this query
  my $key = uc(join($slice->name, $constraint, $logic_name));
  return $self->{'_slice_feature_cache'}{$key} 
    if $self->{'_slice_feature_cache'}{$key};
    
  my $slice_start  = $slice->chr_start();
  my $slice_end    = $slice->chr_end();
  my $slice_strand = $slice->strand();
		 
  my $mapper = 
    $self->db->get_AssemblyMapperAdaptor->fetch_by_type($slice->assembly_type);

  #get the list of contigs this slice is on
  my @cids = 
    $mapper->list_contig_ids( $slice->chr_name, $slice_start ,$slice_end );
  
  return [] unless scalar(@cids);

  my $cid_list = join(',',@cids);

  #get the synonym of the primary_table
  my @tabs = $self->_tables;
  my $syn = $tabs[0]->[1];

  #construct the SQL constraint for the contig ids 
  if($constraint) {
    $constraint .= " AND ${syn}.contig_id IN ($cid_list)";
  } else {
    $constraint = "${syn}.contig_id IN ($cid_list)";
  }

  #for speed the remapping to slice may be done at the time of object creation
  my $features = 
    $self->generic_fetch($constraint, $logic_name, $mapper, $slice); 
  
  if(@$features && (!$features->[0]->can('contig') || 
		    $features->[0]->contig == $slice)) {
    #features have been converted to slice coords already, cache and return
    return $self->{'_slice_feature_cache'}{$key} = $features;
  }

  #remapping has not been done, we have to do our own conversion from
  # raw contig coords to slice coords

  my @out = ();
  
  my ($feat_start, $feat_end, $feat_strand); 

  foreach my $f (@$features) {
    #since feats were obtained in contig coords, attached seq is a contig
    my $contig_id = $f->contig->dbID();

    my ($chr_name, $start, $end, $strand) = 
      $mapper->fast_to_assembly($contig_id, $f->start(), 
				$f->end(),$f->strand(),"rawcontig");

    # undefined start means gap
    next unless defined $start;     

    # maps to region outside desired area 
    next if ($start > $slice_end) || ($end < $slice_start);  
    
    #shift the feature start, end and strand in one call
    if($slice_strand == -1) {
      $f->move( $slice_end - $end + 1, $slice_end - $start + 1, $strand * -1 );
    } else {
      $f->move( $start - $slice_start + 1, $end - $slice_start + 1, $strand );
    }
    
    $f->contig($slice);
    
    push @out,$f;
  }
  
  #update the cache
  return $self->{'_slice_feature_cache'}{$key} = \@out;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::SeqFeature
  Example    : $adaptor->store(@feats);
  Description: ABSTRACT  Subclasses are responsible for implementing this 
               method.  It should take a list of features and store them in 
               the database.
  Returntype : none
  Exceptions : thrown method is not implemented by subclass
  Caller     : general

=cut

sub store{
  my $self = @_;

  $self->throw("Abstract method store not defined by implementing subclass\n");
}


=head2 remove

  Arg [1]    : A feature $feature 
  Example    : $feature_adaptor->remove($feature);
  Description: This removes a feature from the database.  The table the
               feature is removed from is defined by the abstract method
               _tablename, and the primary key of the table is assumed
               to be _tablename() . '_id'.  The feature argument must 
               be an object implementing the dbID method, and for the
               feature to be removed from the database a dbID value must
               be returned.
  Returntype : none
  Exceptions : thrown if $feature arg does not implement dbID(), or if 
               $feature->dbID is not a true value               
  Caller     : general

=cut


sub remove {
  my ($self, $feature) = @_;

  unless($feature->can('dbID')) {
    $self->throw("Feature [$feature] does not implement method dbID");
  }

  unless($feature->dbID) {
    $self->warn("BaseFeatureAdaptor::remove - dbID not defined - " .
                "feature could not be removed");
  }

  my @tabs = $self->_tables;
  my ($table) = @{$tabs[0]};

  my $sth = $self->prepare("DELETE FROM $table WHERE ${table}_id = ?");
  $sth->execute($feature->dbID());

  #unset the feature dbID
  $feature->dbID('');
  
  return;
}



=head2 remove_by_RawContig

  Arg [1]    : Bio::EnsEMBL::RawContig $contig 
  Example    : $feature_adaptor->remove_by_RawContig($contig);
  Description: This removes features from the database which lie on a removed
               contig.  The table the features are removed from is defined by 
               the abstract method_tablename, and the primary key of the table
               is assumed to be contig_id.
  Returntype : none
  Exceptions : thrown if no contig is supplied
  Caller     : general

=cut

sub remove_by_RawContig {
  my ($self, $contig) = @_;

  unless($contig) {
    $self->throw("BaseFeatureAdaptor::remove - no contig supplied: ".
		 "Deletion of features failed.");
  }

  my @tabs = $self->_tables;

  my ($table_name) = @{$tabs[0]};

  my $sth = $self->prepare("DELETE FROM $table_name
                            WHERE contig_id = ?");

  $sth->execute($contig->dbID);

  return;
}



=head2 _tables

  Args       : none
  Example    : $tablename = $self->_table_name()
  Description: ABSTRACT PROTECTED Subclasses are responsible for implementing
               this method.  It should list of [tablename, alias] pairs.  
               Additionally the primary table (with the dbID, analysis_id, and
               score) should be the first table in the list.
               e.g:
               ( ['repeat_feature',   'rf'],
                 ['repeat_consensus', 'rc']);
               used to obtain features.  
  Returntype : list of [tablename, alias] pairs
  Exceptions : thrown if not implemented by subclass
  Caller     : BaseFeatureAdaptor::generic_fetch

=cut

sub _tables {
  my $self = shift;

  $self->throw("abstract method _tables not defined by implementing" .
               " subclass of BaseFeatureAdaptor");
  return undef;
}


=head2 _columns

  Args       : none
  Example    : $tablename = $self->_columns()
  Description: ABSTRACT PROTECTED Subclasses are responsible for implementing
               this method.  It should return a list of columns to be used
               for feature creation
  Returntype : list of strings
  Exceptions : thrown if not implemented by subclass
  Caller     : BaseFeatureAdaptor::generic_fetch

=cut

sub _columns {
  my $self = shift;

  $self->throw("abstract method _columns not defined by implementing" .
               " subclass of BaseFeatureAdaptor");
}



=head2 _left_join

  Arg [1]    : none
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut





=head2 _default_where_clause

  Arg [1]    : none
  Example    : none
  Description: May be overridden to provide an additional where constraint to 
               the SQL query which is generated to fetch feature records.
               This constraint is always appended to the end of the generated
               where clause
  Returntype : string
  Exceptions : none
  Caller     : generic_fetch

=cut

sub _default_where_clause {
  my $self = shift;

  return '';
}



=head2 _left_join

  Arg [1]    : none
  Example    : none
  Description: Can be overridden by a subclass to specify any left joins
               which should occur. The table name specigfied in the join
               must still be present in the return values of 
  Returntype : a {'tablename' => 'join condition'} pair 
  Exceptions : none
  Caller     : general

=cut

sub _left_join {
  my $self = shift;

  return '';
}



=head2 _final_clause

  Arg [1]    : none
  Example    : none
  Description: May be overriden to provide an additional clause to the end
               of the SQL query used to fetch feature records.  
               This is useful to add a required ORDER BY clause to the 
               query for example.
  Returntype : string
  Exceptions : none
  Caller     : generic_fetch

=cut

sub _final_clause {
  my $self = shift;

  return '';
}

=head2 _objs_from_sth

  Arg [1]    : DBI::row_hashref $hashref containing key-value pairs 
               for each of the columns specified by the _columns method
  Example    : my @feats = $self->_obj_from_hashref
  Description: ABSTRACT PROTECTED The subclass is responsible for implementing
               this method.  It should take in a DBI row hash reference and
               return a list of created features in contig coordinates.
  Returntype : list of Bio::EnsEMBL::*Features in contig coordinates
  Exceptions : thrown if not implemented by subclass
  Caller     : BaseFeatureAdaptor::generic_fetch

=cut

sub _objs_from_sth {
  my $self = shift;

  $self->throw("abstract method _obj_from_sth not defined by implementing"
             . " subclass of BaseFeatureAdaptor");
} 


=head2 deleteObj

  Arg [1]    : none
  Example    : none
  Description: Cleans up internal caches and references to other objects so
               that correct garbage collection may occur.
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBConnection::deleteObj

=cut

sub deleteObj {
  my $self = shift;

  #flush feature cache
  %{$self->{'_slice_feature_cache'}} = ();
}

1;



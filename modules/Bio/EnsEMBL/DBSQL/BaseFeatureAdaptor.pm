#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor
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

Abstract class should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for feature adaptors. This base class is simply a way
of eliminating code duplication through the implementation of methods 
common to all feature adaptors.

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Cache;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

my $SLICE_FEATURE_CACHE_SIZE = 12;


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
  
  my $tablename = $self->_tablename();
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
    
    if($constraint) {
      $constraint .= " AND analysis_id = $analysis_id";
    } else {
      $constraint = " analysis_id = $analysis_id";
    }
  } 
      
  my $sql = "SELECT $columns FROM $tablename " . 
    ($constraint ? " where $constraint " : '' );
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

  my $tablename = $self->_tablename();
  my $constraint = "${tablename}_id = $id";

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

  if($constraint) {
    $constraint .= " AND contig_id = $cid";
  } else {
    $constraint = "contig_id = $cid";
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
    $constraint = "score > $score";
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
    $constraint = "score > $score";
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

  #check the cache and return if we have already done this query
  $logic_name = '' unless defined $logic_name;

  my $key = join($slice->name, $constraint, $logic_name);
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

  #construct the SQL constraint for the contig ids 
  if($constraint) {
    $constraint .= " AND contig_id IN ($cid_list)";
  } else {
    $constraint = "contig_id IN ($cid_list)";
  }

  #for speed the remapping to slice may be done at the time of object creation
  my $features = 
    $self->generic_fetch($constraint, $logic_name, $mapper, $slice); 
  
  if(@$features && $features->[0]->contig == $slice) {
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
               feature to be removed from the datasbase a dbID value must
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

  my $table = $self->_tablename();

  my $sth = $self->prepare("DELETE FROM $table WHERE ${table}_id = ?");
  $sth->execute($feature->dbID());

  #unset the feature dbID
  $feature->dbID('');
  
  return;
}


=head2 _tablename

  Args       : none
  Example    : $tablename = $self->_table_name()
  Description: ABSTRACT PROTECTED Subclasses are responsible for implementing
               this method.  It should return the name of the table to be
               used to obtain features.  
  Returntype : string
  Exceptions : thrown if not implemented by subclass
  Caller     : BaseFeatureAdaptor::generic_fetch

=cut

sub _tablename {
  my $self = shift;

  $self->throw("abstract method _tablename not defined by implementing" .
               " subclass of AlignFeatureAdaptor");
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
               " subclass of AlignFeatureAdaptor");
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

  $self->throw("abstract method _obj_from_hashref not defined by implementing"
             . " subclass of AlignFeatureAdaptor");
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




=head2 fetch_by_Contig_constraint

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_RawContig_constraint instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Contig_constraint {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Contig_constraint has been renamed fetch_all_by_RawContig_constraint\n" . caller);

  return $self->fetch_all_by_RawContig_constraint(@args);
}



=head2 fetch_by_Contig

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_RawContig instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Contig {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Contig has been renamed fetch_all_by_RawContig\n" . caller);

  return $self->fetch_all_by_RawContig(@args);
}


=head2 fetch_by_Contig_and_score

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_RawContig_and_score instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Contig_and_score {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Contig_and_score has been renamed fetch_all_by_RawContig_and_score\n" . caller);

  return $self->fetch_all_by_RawContig_and_score(@args);
}



=head2 fetch_all_by_Contig

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_RawContig instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_all_by_Contig {
  my ($self, @args) = @_;

  $self->warn("fetch_all_by_Contig has been renamed fetch_all_by_RawContig\n" . caller);

  return $self->fetch_all_by_RawContig(@args);
}


=head2 fetch_all_by_Contig_and_score

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_RawContig_and_score instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_all_by_Contig_and_score {
  my ($self, @args) = @_;

  $self->warn("fetch_all_by_Contig_and_score has been renamed fetch_all_by_RawContig_and_score\n" . caller);

  return $self->fetch_all_by_RawContig_and_score(@args);
}


=head2 fetch_all_by_Contig_constraint

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_RawContig_constraint instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_all_by_Contig_constraint {
  my ($self, @args) = @_;

  $self->warn("fetch_all_by_Contig_constraint has been renamed fetch_all_by_RawContig_constraint\n" . caller);

  return $self->fetch_all_by_RawContig_constraint(@args);
}




=head2 fetch_by_Slice_and_score

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_Slice_and_score instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Slice_and_score {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Slice_and_score has been renamed fetch_all_by_Slice_and_score\n" . caller);

  return $self->fetch_all_by_Slice_and_score(@args);
}


=head2 fetch_by_Slice_constraint

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_Slice_constraint instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Slice_constraint {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Slice_constraint has been renamed fetch_all_by_Slice_constraint\n" . caller);

  return $self->fetch_all_by_Slice_constraint(@args);
}



=head2 fetch_by_Slice

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_Slice instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Slice {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Slice has been renamed fetch_all_by_Slice\n" . caller);

  return $self->fetch_all_by_Slice(@args);
}




1;



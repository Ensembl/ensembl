=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::BaseAdaptor - Base Adaptor for DBSQL adaptors

=head1 SYNOPSIS

  # base adaptor provides

  # SQL prepare function
  $adaptor->prepare("sql statement");

  # get of root DBAdaptor object
  $adaptor->db();

  # constructor, ok for inheritence
  $adaptor = Bio::EnsEMBL::DBSQL::SubClassOfBaseAdaptor->new($dbobj)

=head1 DESCRIPTION

This is a true base class for Adaptors in the Ensembl DBSQL
system. Original idea from Arne

Adaptors are expected to have the following functions

  $obj = $adaptor->fetch_by_dbID($internal_id);

which builds the object from the primary key of the object. This
function is crucial because it allows adaptors to collaborate relatively
independently of each other - in other words, we can change the schema
under one adaptor without too many knock on changes through the other
adaptors.

Most adaptors will also have

  $dbid = $adaptor->store($obj);

which stores the object. Currently the storing of an object also causes
the objects to set

  $obj->dbID();

correctly and attach the adaptor.

Other fetch functions go by the convention of

  @object_array = @{ $adaptor->fetch_all_by_XXXX($arguments_for_XXXX) };

sometimes it returns an array ref denoted by the 'all' in the name of
the method, sometimes an individual object. For example

  $gene = $gene_adaptor->fetch_by_stable_id($stable_id);

or

  @fp = @{ $simple_feature_adaptor->fetch_all_by_Slice($slice) };

Occassionally adaptors need to provide access to lists of ids. In this
case the convention is to go list_XXXX, such as

  @gene_ids = @{ $gene_adaptor->list_geneIds() };

(note: this method is poorly named)

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::BaseAdaptor;
require Exporter;
use vars qw(@ISA @EXPORT);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use DBI qw(:sql_types);
use Data::Dumper;

@ISA = qw(Exporter);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection $dbobj
  Example    : $adaptor = new AdaptorInheritedFromBaseAdaptor($dbobj);
  Description: Creates a new BaseAdaptor object.  The intent is that this
               constructor would be called by an inherited superclass either
               automatically or through $self->SUPER::new in an overridden 
               new method.
  Returntype : Bio::EnsEMBL::DBSQL::BaseAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBConnection
  Status     : Stable

=cut

sub new {
  my ( $class, $dbobj ) = @_;

  my $self = bless {}, $class;

  if ( !defined $dbobj || !ref $dbobj ) {
    throw("Don't have a db [$dbobj] for new adaptor");
  }

  if ( $dbobj->isa('Bio::EnsEMBL::DBSQL::DBAdaptor') ) {
    $self->db($dbobj);
    $self->dbc( $dbobj->dbc );
    $self->species_id( $dbobj->species_id() );
    $self->is_multispecies( $dbobj->is_multispecies() );
  } elsif ( ref($dbobj) =~ /DBAdaptor$/ ) {
    $self->db($dbobj);
    $self->dbc( $dbobj->dbc );
  } elsif ( ref($dbobj) =~ /DBConnection$/ ) {
    $self->dbc($dbobj);
  } else {
    throw("Don't have a DBAdaptor [$dbobj] for new adaptor");
  }

  return $self;
}


=head2 prepare

  Arg [1]    : string $string
               a SQL query to be prepared by this adaptors database
  Example    : $sth = $adaptor->prepare("select yadda from blabla")
  Description: provides a DBI statement handle from the adaptor. A convenience
               function so you dont have to write $adaptor->db->prepare all the
               time
  Returntype : DBI::StatementHandle
  Exceptions : none
  Caller     : Adaptors inherited from BaseAdaptor
  Status     : Stable

=cut

sub prepare {
  my ( $self, $string ) = @_;

  # Uncomment next line to cancel caching on the SQL side.
  # Needed for timing comparisons etc.
  #$string =~ s/SELECT/SELECT SQL_NO_CACHE/i;

  return $self->dbc->prepare($string);
}


=head2 db

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor $obj 
               the database this adaptor is using.
  Example    : $db = $adaptor->db();
  Description: Getter/Setter for the DatabaseConnection that this adaptor is 
               using.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : Adaptors inherited from BaseAdaptor
  Status     : Stable

=cut

sub db {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'db'} = $value;
  }

  return $self->{'db'};
}

=head2 dbc

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::DBConnection $obj 
               the database this adaptor is using.
  Example    : $db = $adaptor->db();
  Description: Getter/Setter for the DatabaseConnection that this adaptor is 
               using.
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : none
  Caller     : Adaptors inherited from BaseAdaptor
  Status     : Stable

=cut

sub dbc {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'dbc'} = $value;
  }

  return $self->{'dbc'};
}

=head2 is_multispecies

  Arg [1]    : (optional) boolean $arg
  Example    : if ($adaptor->is_multispecies()) { }
  Description: Getter/Setter for the is_multispecies boolean of
               to use for this adaptor.
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_multispecies {
  my ( $self, $arg ) = @_;

  if ( defined($arg) ) {
    $self->{_is_multispecies} = $arg;
  }

  return $self->{_is_multispecies};
}

=head2 species_id

  Arg [1]    : (optional) int $species_id
               The internal ID of the species in a multi-species database.
  Example    : $db = $adaptor->db();
  Description: Getter/Setter for the internal ID of the species in a
               multi-species database.  The default species ID is 1.
  Returntype : Integer
  Exceptions : none
  Caller     : Adaptors inherited from BaseAdaptor
  Status     : Stable

=cut

sub species_id {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'species_id'} = $value;
  }

  return $self->{'species_id'} || 1;
}


# list primary keys for a particular table
# args are table name and primary key field
# if primary key field is not supplied, tablename_id is assumed
# returns listref of IDs
sub _list_dbIDs {
  my ( $self, $table, $pk, $ordered ) = @_;

  if ( !defined($pk) ) { $pk = $table . "_id" }

  my $sql = "SELECT " . $pk . "  FROM " . $table;

  if ( defined($ordered) && $ordered ) {
    $sql .= " order by seq_region_id, seq_region_start";
  }

  my $sth = $self->prepare($sql);
  eval { $sth->execute() };
  if ($@) {
    throw("Detected an error whilst executing SQL '${sql}': $@");
  }

  my @out;
  while ( my ($id) = $sth->fetchrow() ) { push( @out, $id ) }

  $sth->finish();

  return \@out;
}


# _straight_join

#   Arg [1]    : (optional) boolean $new_val
#   Example    : $self->_straight_join(1);
#                $self->generic_fetch($constraint);
#                $self->_straight_join(0);
#   Description: PROTECTED Getter/Setter that turns on/off the use of 
#                a straight join in queries.
#   Returntype : boolean
#   Exceptions : none
#   Caller     : general

sub _straight_join {
  my $self = shift;
  if(@_) {
    $self->{'_straight_join'} = shift;
  }

  return $self->{'_straight_join'};
}


=head2 bind_param_generic_fetch

 Arg [1]   : (optional)  scalar $param
              This is the parameter to bind
 Arg [2]   : (optional) int $sql_type
              Type of the parameter (from DBI (:sql_types))
 Example   :  $adaptor->bind_param_generic_fetch($stable_id,SQL_VARCHAR);
              $adaptor->generic_fetch();
 Description:  When using parameters for the query, will call the bind_param to avoid
               some security issues. If there are no arguments, will return the bind_parameters
 ReturnType : listref
 Exceptions:  if called with one argument

=cut

sub bind_param_generic_fetch{
    my $self = shift;
    my $param = shift;
    my $sql_type = shift;

    if (defined $param && !defined $sql_type){
	throw("Need to specify sql_type for parameter $param\n");
    }
    elsif (defined $param && defined $sql_type){
	#check when there is a SQL_INTEGER type that the parameter is really a number
	if ($sql_type eq SQL_INTEGER){
	    throw "Trying to assign a non numerical parameter to an integer value in the database" if ($param !~ /^\d+$/);
	}
	#both paramters have been entered, push it to the bind_param array
	push @{$self->{'_bind_param_generic_fetch'}},[$param,$sql_type];
    }
    elsif (!defined $param && !defined $sql_type){
	#when there are no arguments, return the array
	return $self->{'_bind_param_generic_fetch'};
    }
	
}



=head2 generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) Bio::EnsEMBL::AssemblyMapper $mapper
               A mapper object used to remap features
               as they are retrieved from the database
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice that features should be remapped to
  Example    : $fts = $a->generic_fetch('contig_id in (1234, 1235)', 'Swall');
  Description: Performs a database fetch and returns feature objects in
               contig coordinates.
  Returntype : listref of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : none
  Caller     : BaseFeatureAdaptor, ProxyDnaAlignFeatureAdaptor::generic_fetch
  Status     : Stable

=cut

sub generic_fetch {
  my ($self, $constraint, $mapper, $slice) = @_;

  my @tabs = $self->_tables();

  my $extra_default_where;

  # Hack for feature types that needs to be restricted to species_id (in
  # coord_system).
  if (    $self->is_multispecies()
       && $self->isa('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor')
       && !$self->isa('Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor') )
  {
    # We do a check to see if there is already seq_region
    # and coord_system defined to ensure we get the right
    # alias.  We then do the extra query irrespectively of
    # what has already been specified by the user.
    my %thash = map { $_->[0] => $_->[1] } @tabs;

    my $sr_alias =
      ( exists( $thash{seq_region} ) ? $thash{seq_region} : 'sr' );
    my $cs_alias =
      ( exists( $thash{coord_system} ) ? $thash{coord_system} : 'cs' );

    if ( !exists( $thash{seq_region} ) ) {
      push( @tabs, [ 'seq_region', $sr_alias ] );
    }
    if ( !exists( $thash{coord_system} ) ) {
      push( @tabs, [ 'coord_system', $cs_alias ] );
    }

    $extra_default_where = sprintf(
                      '%s.seq_region_id = %s.seq_region_id '
                        . 'AND %s.coord_system_id = %s.coord_system_id '
                        . 'AND %s.species_id = ?',
                      $tabs[0]->[1], $sr_alias, $sr_alias,
                      $cs_alias,     $cs_alias );

    $self->bind_param_generic_fetch( $self->species_id(), SQL_INTEGER );
  } ## end if ( $self->is_multispecies...)

  my $columns = join(', ', $self->_columns());

  my $db = $self->db();

  #
  # Construct a left join statement if one was defined, and remove the
  # left-joined table from the table list
  #
  my @left_join_list = $self->_left_join();
  my $left_join_prefix = '';
  my $left_join = '';
  my @tables;
  if(@left_join_list) {
    my %left_join_hash = map { $_->[0] => $_->[1] } @left_join_list;
    while(my $t = shift @tabs) {
      if( exists $left_join_hash{ $t->[0] } ) {
        my $condition = $left_join_hash{ $t->[0] };
        my $syn = $t->[1];
        $left_join .=
          "\n  LEFT JOIN " . $t->[0] . " $syn ON $condition ) ";
        $left_join_prefix .= '(';
      } else {
        push @tables, $t;
      }
    }
  } else {
    @tables = @tabs;
  }

  my $straight_join = '';

  if($self->_straight_join()) {
    $straight_join = "STRAIGHT_JOIN";
  }

  #construct a nice table string like 'table1 t1, table2 t2'
  my $tablenames = join(', ', map({ join(' ', @$_) } @tables));

  my $sql =
      "SELECT $straight_join $columns\n"
    . "FROM $left_join_prefix ($tablenames) $left_join";

  my $default_where = $self->_default_where_clause();
  my $final_clause = $self->_final_clause;

  if ($extra_default_where) {
    if ($default_where) {
      $default_where .= "\n AND $extra_default_where";
    } else {
      $default_where = $extra_default_where;
    }
  }

  #append a where clause if it was defined
  if ($constraint) {
    $sql .= "\n WHERE $constraint ";
    if ($default_where) {
      $sql .= " AND\n       $default_where ";
    }
  } elsif ($default_where) {
    $sql .= "\n WHERE $default_where ";
  }

  #append additional clauses which may have been defined
  $sql .= "\n$final_clause";


  # FOR DEBUG:
  # printf(STDERR "SQL:\n%s\n", $sql);


  my $sth = $db->dbc->prepare($sql);
  my $bind_parameters = $self->bind_param_generic_fetch();
  if (defined $bind_parameters){
      #if we have bind the parameters, call the DBI to bind them
      my $i = 1;
      foreach my $param (@{$bind_parameters}){
	  $sth->bind_param($i,$param->[0],$param->[1]);
	  $i++;
      }
      #after binding parameters, undef for future queries
      $self->{'_bind_param_generic_fetch'} = ();
  }
  eval { $sth->execute() };
  if ($@) {
    throw("Detected an error whilst executing SQL '${sql}': $@");
  }

  my $res = $self->_objs_from_sth($sth, $mapper, $slice);
  $sth->finish();
  return $res;
}


=head2 fetch_by_dbID

  Arg [1]    : int $id
               The unique database identifier for the feature to be obtained
  Example    : $feat = $adaptor->fetch_by_dbID(1234));
               $feat = $feat->transform('contig');
  Description: Returns the feature created from the database defined by the
               the id $id.  The feature will be returned in its native
               coordinate system.  That is, the coordinate system in which it
               is stored in the database.  In order to convert it to a
               particular coordinate system use the transfer() or transform()
               method.  If the feature is not found in the database then
               undef is returned instead
  Returntype : Bio::EnsEMBL::Feature or undef
  Exceptions : thrown if $id arg is not provided
               does not exist
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID{
  my ($self,$id) = @_;

  throw("id argument is required") if(!defined $id);

  #construct a constraint like 't1.table1_id = 123'
  my @tabs = $self->_tables;
  my ($name, $syn) = @{$tabs[0]};
  $self->bind_param_generic_fetch($id,SQL_INTEGER);
  my $constraint = "${syn}.${name}_id = ?";

  #Should only be one
  my ($feat) = @{$self->generic_fetch($constraint)};

  return undef if(!$feat);

  return $feat;
}


=head2 fetch_all_by_dbID_list

  Arg [1]    : listref of integers $id_list
               The unique database identifiers for the features to
               be obtained.
  Example    : @feats = @{$adaptor->fetch_all_by_dbID_list([1234, 2131, 982]))};
  Description: Returns the features created from the database
               defined by the the IDs in contained in the provided
               ID list $id_list.  The features will be returned
               in their native coordinate system.  That is, the
               coordinate system in which they are stored in the
               database.  In order to convert the features to a
               particular coordinate system use the transfer() or
               transform() method.  If none of the features are
               found in the database a reference to an empty list is
               returned.
  Returntype : listref of Bio::EnsEMBL::Features
  Exceptions : thrown if $id arg is not provided
               does not exist
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_dbID_list {
  my ( $self, $id_list_ref ) = @_;

  if ( !defined($id_list_ref) || ref($id_list_ref) ne 'ARRAY' ) {
    throw("id_list list reference argument is required");
  }

  if ( !@{$id_list_ref} ) { return [] }

  # Construct a constraint like 't1.table1_id = 123'
  my @tabs = $self->_tables();
  my ( $name, $syn ) = @{ $tabs[0] };

  # Ensure that we do not exceed MySQL's max_allowed_packet (defaults to
  # 1 MB) splitting large queries into smaller queries of at most 256 KB
  # (32768 8-bit characters).  Assuming a (generous) average dbID string
  # length of 16, this means 2048 dbIDs in each query.
  my $max_size = 2048;

  my @id_list = @{$id_list_ref};

  my @out;

  while (@id_list) {
    my @ids;
    my $id_str;

    if ( scalar(@id_list) > $max_size ) {
      @ids = splice( @id_list, 0, $max_size );
    } else {
      @ids     = @id_list;
      @id_list = ();
    }

    if ( scalar(@ids) > 1 ) {
      $id_str = " IN (" . join( ',', @ids ) . ")";
    } else {
      $id_str = " = " . $ids[0];
    }

    my $constraint = "${syn}.${name}_id $id_str";

    push @out, @{ $self->generic_fetch($constraint) };
  }

  return \@out;
} ## end sub fetch_all_by_dbID_list

# might not be a good idea, but for convenience
# shouldnt be called on the BIG tables though

sub fetch_all {
  my $self = shift;
  return $self->generic_fetch();
}


#_tables
#
#  Args       : none
#  Example    : $tablename = $self->_table_name()
#  Description: ABSTRACT PROTECTED
#               Subclasses are responsible for implementing this
#               method.  It should list of [tablename, alias] pairs.
#               Additionally the primary table (with the dbID,
#               analysis_id, and score) should be the first table in
#               the list. e.g:
#               ( ['repeat_feature',   'rf'],
#                 ['repeat_consensus', 'rc']);
#               used to obtain features.  
#  Returntype : list of [tablename, alias] pairs
#  Exceptions : thrown if not implemented by subclass
#  Caller     : BaseFeatureAdaptor::generic_fetch
#

sub _tables {
  throw(   "abstract method _tables not defined "
         . "by implementing subclass of BaseAdaptor" );
}


#_columns
#
#  Args       : none
#  Example    : $tablename = $self->_columns()
#  Description: ABSTRACT PROTECTED
#               Subclasses are responsible for implementing this
#               method.  It should return a list of columns to be
#               used for feature creation.
#  Returntype : list of strings
#  Exceptions : thrown if not implemented by subclass
#  Caller     : BaseFeatureAdaptor::generic_fetch
#

sub _columns {
  throw(   "abstract method _columns not defined "
         . "by implementing subclass of BaseAdaptor" );
}


# _default_where_clause
#
#  Arg [1]    : none
#  Example    : none
#  Description: May be overridden to provide an additional where
#               constraint to the SQL query which is generated to
#               fetch feature records.  This constraint is always
#               appended to the end of the generated where clause
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch
#

sub _default_where_clause { return '' }


# _left_join

#  Arg [1]    : none
#  Example    : none
#  Description: Can be overridden by a subclass to specify any left
#               joins which should occur.  The table name specigfied
#               in the join must still be present in the return
#               values of.
#  Returntype : a {'tablename' => 'join condition'} pair
#  Exceptions : none
#  Caller     : general
#

sub _left_join { return () }


#_final_clause

#  Arg [1]    : none
#  Example    : none
#  Description: May be overriden to provide an additional clause
#               to the end of the SQL query used to fetch feature
#               records.  This is useful to add a required ORDER BY
#               clause to the query for example.
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch

sub _final_clause { return '' }


#_objs_from_sth

#  Arg [1]    : DBI::row_hashref $hashref containing key-value pairs 
#               for each of the columns specified by the _columns method
#  Example    : my @feats = $self->_obj_from_hashref
#  Description: ABSTRACT PROTECTED
#               The subclass is responsible for implementing this
#               method.  It should take in a DBI row hash reference
#               and return a list of created features in contig
#               coordinates.
#  Returntype : list of Bio::EnsEMBL::*Features in contig coordinates
#  Exceptions : thrown if not implemented by subclass
#  Caller     : BaseFeatureAdaptor::generic_fetch

sub _objs_from_sth {
  throw(   "abstract method _objs_from_sth not defined "
         . "by implementing subclass of BaseAdaptor" );
}

sub dump_data {
  my $self = shift;
  my $data = shift;

  my $dumper = Data::Dumper->new([$data]);
  $dumper->Indent(0);
  $dumper->Terse(1);
   my $dump = $dumper->Dump();
# $dump =~ s/'/\\'/g; 
 # $dump =~ s/^\$VAR1 = //;
  return $dump;
}

sub get_dumped_data {
    my $self = shift;
    my $data = shift;

    $data =~ s/\n|\r|\f|\\//g;
    return eval ($data);
}


1;

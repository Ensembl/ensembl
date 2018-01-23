=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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
  $adaptor = Bio::EnsEMBL::DBSQL::SubClassOfBaseAdaptor->new($dbobj);

=head1 DESCRIPTION

This is a true base class for Adaptors in the Ensembl DBSQL
system.

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

=cut

package Bio::EnsEMBL::DBSQL::BaseAdaptor;
require Exporter;
use vars qw(@ISA @EXPORT);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref assert_integer wrap_array);
use DBI qw(:sql_types);
use Data::Dumper;
use Scalar::Util qw/looks_like_number/;

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

  my $sql = sprintf( "SELECT `%s` FROM `%s`", $pk, $table );

  my $join_with_cs = 0;
  if (    $self->is_multispecies()
       && $self->isa('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor')
       && !$self->isa('Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor') )
  {
    
    $sql .= q(
JOIN seq_region USING (seq_region_id)
JOIN coord_system cs USING (coord_system_id)
WHERE cs.species_id = ?
);

    $join_with_cs = 1;
  }

  if ( defined($ordered) && $ordered ) {
    $sql .= " ORDER BY seq_region_id, seq_region_start";
  }

  my $sth = $self->prepare($sql);

  if ($join_with_cs) {
    $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
  }

  eval { $sth->execute() };
  if ($@) {
    throw("Detected an error whilst executing SQL '${sql}': $@");
  }

  my $id;
  $sth->bind_col( 1, \$id );

  my @out;
  while ( $sth->fetch() ) {
    push( @out, $id );
  }

  return \@out;
} ## end sub _list_dbIDs


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

sub _can_straight_join {
    my $self = shift;
    return $self->dbc->_driver_object->can_straight_join;
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
	    throw "Trying to assign a non numerical parameter to an integer value in the database" if ($param !~ /^[+-]{0,1}\d+$/);
	}
	#both paramters have been entered, push it to the bind_param array
	push @{$self->{'_bind_param_generic_fetch'}},[$param,$sql_type];
    }
    elsif (!defined $param && !defined $sql_type){
	#when there are no arguments, return the array
	return $self->{'_bind_param_generic_fetch'};
    }
	
}

# Used to reset the params without circumventing scope
sub _bind_param_generic_fetch {
  my ($self, $_bind_param_generic_fetch) = @_;
  $self->{'_bind_param_generic_fetch'} = $_bind_param_generic_fetch if $_bind_param_generic_fetch;
  return $self->{_bind_param_generic_fetch};
}

=head2 generate_in_constraint
  
  Arg [1]     : ArrayRef or Scalar $list
                List or a single value of items to be pushed into an IN statement
  Arg [2]     : Scalar $column
                Column this IN statement is being applied to. Please fully resolve the
                column.
  Arg [3]     : Scalar $param_type
                Data type which should be used when binding. Please use DBI data type symbols
  Arg [4]     : Scalar boolean $inline_variables
                Boolean to control if variables are inlined in the constraint. If
                false values are bound via bind_param_generic_fetch() (the default behaviour).

  Description : Used internally to generate a SQL constraint to restrict a query by an IN statement.
                The code generates the complete IN statement.
  Returntype  : String
  Exceptions  : If no list is supplied, the list of values is empty or no data type was given
  Caller      : general

=cut

sub generate_in_constraint {
  my ($self, $list, $column, $param_type, $inline_variables) = @_;
  throw("A list of values must be given") if ! defined $list;
  $list = wrap_array($list); # homogenise into an array
  throw "We should be given at least one value to insert" if scalar(@{$list}) == 0;
  throw "Please supply the DBI param type" if ! defined $param_type;
  #Figure out if we need to quote our values if we are asked to inline the variables
  my $quote_values = 1;
  if($param_type == SQL_INTEGER || $param_type == SQL_TINYINT || $param_type == SQL_DOUBLE ) {
    $quote_values = 0;
  }

  my $constraint = qq{${column} IN (};
  if($inline_variables) {
    if($quote_values) {
      $constraint .= join(q{,}, map { qq{"${_}"} } @{$list});  
    }
    else {
      $constraint .= join(q{,}, @{$list});  
    }
  }
  else {
    my @subs = ('?') x scalar(@{$list});
    $constraint .= join(q{,}, @subs);
    $self->bind_param_generic_fetch($_, $param_type) for @{$list};
  }
  $constraint .= q{)};
  return $constraint;
}

=head2 generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) Bio::EnsEMBL::AssemblyMapper $mapper
               A mapper object used to remap features
               as they are retrieved from the database
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice that features should be remapped to
  Example    : $fts = $a->generic_fetch('contig_id in (1234, 1235)');
  Description: Performs a database fetch and returns feature objects in
               contig coordinates.
  Returntype : listref of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : Thrown if there is an issue with querying the data
  Caller     : BaseFeatureAdaptor, ProxyDnaAlignFeatureAdaptor::generic_fetch
  Status     : Stable

=cut

sub generic_fetch {
  my ($self, $constraint, $mapper, $slice) = @_;
  my $sql = $self->_generate_sql($constraint);
  my $params = $self->bind_param_generic_fetch();
  $params ||= [];
  $self->{_bind_param_generic_fetch} = undef;
  my $sth = $self->db()->dbc()->prepare($sql);
  my $i = 1;
  foreach my $param (@{$params}){
    $sth->bind_param($i,$param->[0],$param->[1]);
    $i++;
  }
  eval { $sth->execute() };
  if ($@) {
    throw("Detected an error whilst executing SQL '${sql}': $@");
  }

  my $res = $self->_objs_from_sth($sth, $mapper, $slice);
  $sth->finish();
  return $res;
}

=head2 generic_count

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Example    : $number_feats = $a->generic_count('contig_id in (1234, 1235)');
  Description: Performs a database fetch and returns a count of those features
               found. This is analagous to C<generic_fetch()>
  Returntype : Integer count of the elements.
  Exceptions : Thrown if there is an issue with querying the data

=cut

sub generic_count {
  my ($self, $constraint) = @_;
  my $sql = $self->_generate_sql($constraint, 'count(*)');
  my $params = $self->bind_param_generic_fetch();
  $params ||= [];
  $self->{_bind_param_generic_fetch} = undef;
  my $h = $self->db()->dbc()->sql_helper();
  my $count = $h->execute_single_result(-SQL => $sql, -PARAMS => $params);
  return $count;
}

sub _generate_sql {
  my ($self, $constraint, @input_columns) = @_;
  
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

  @input_columns = $self->_columns() if ! @input_columns;
  my $columns = join(', ', @input_columns);

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
        my $t_alias = $t->[0] . " " . $t->[1];
      if( exists $left_join_hash{ $t->[0] } || exists $left_join_hash{$t_alias}) {
        my $condition = $left_join_hash{ $t->[0] };
        $condition ||= $left_join_hash{$t_alias};
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

  if($self->_straight_join() and $self->_can_straight_join) {
    $straight_join = "STRAIGHT_JOIN";
  }

  #construct a nice table string like 'table1 t1, table2 t2'
  my $tablenames = join(', ', map({ join(' ', @$_) } @tables));

  my $sql =
      "SELECT $straight_join $columns \n"
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
  #printf(STDERR "SQL:\n%s\n", $sql);
  
  return $sql;
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

sub fetch_by_dbID {
  my ($self, $id) = @_;
  if ($self->_no_id_cache()) {
    return $self->_uncached_fetch_by_dbID($id);
  }
  return $self->_id_cache()->get($id);
}

# The actual implmenetation moved sideways to allow for uncached access
# otherwise we'd constantly loop

sub _uncached_fetch_by_dbID{
  my ($self,$id) = @_;

  throw("id argument is required") if(!defined $id);

  #construct a constraint like 't1.table1_id = 123'
  my @tabs = $self->_tables;
  my ($name, $syn) = @{$tabs[0]};
  $self->bind_param_generic_fetch($id,SQL_INTEGER);
  my $constraint = "${syn}.${name}_id = ?";

  #Should only be one
  my ($feat) = @{$self->generic_fetch($constraint)};

  return if(!$feat);

  return $feat;
}


=head2 fetch_all_by_dbID_list

  Arg [1]    : listref of integers $id_list
               The unique database identifiers for the features to
               be obtained.
  Arg [2]    : optional - Bio::EnsEMBL::Slice to map features onto.
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
  my ($self, $id_list_ref, $slice) = @_;
  if ($self->_no_id_cache()) {
    return $self->_uncached_fetch_all_by_dbID_list($id_list_ref, $slice);
  }
  return $self->_id_cache()->get_by_list($id_list_ref, $slice);
}

# The actual implmenetation moved sideways to allow for uncached access
# otherwise we'd constantly loop
sub _uncached_fetch_all_by_dbID_list {
  my ( $self, $id_list_ref, $slice ) = @_;
  return $self->_uncached_fetch_all_by_id_list($id_list_ref, $slice, "dbID", 1);
} ## end sub fetch_all_by_dbID_list

=head2 _uncached_fetch_all_by_id_list

  Arg [1]    : listref of IDs
  Arg [2]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice that features should be remapped to
  Arg [3]    : String describing the ID type.
               Valid values include dbID and stable_id. dbID is an alias for
               the primary key, while other names map directly to table columns
               of the Feature this adaptor manages.
  Arg [4]    : Boolean $numeric
               Indicates if the incoming data is to be processed as a numeric
               or as a String. If arg [3] was set to dbID then we default this to
               be true. If arg [3] was set to stable_id then we default this to
               be false.
               When not using a standard arg[3] the IDs are assumed to be Strings.
  Arg [5]    : Integer $max_size
               Control the maximum number of IDs sent to a database in a single 
               query. Defaults to 2K for Strings and 16K for integers. Only
               provide if you know *exactly* why you need to edit it.
  Example    : $list_of_features = $adaptor->_uncached_fetch_all_by_id_list(
                   [qw(ENSG00000101321 ENSG00000101346 ENSG00000101367)],
                   undef,
                   "stable_id", 0); #using strings
               
               # Numeric set to true because we are using numerics
               $list_of_features = $adaptor->_uncached_fetch_all_by_id_list(
                   [1,2,3,4],
                   undef,
                   "dbID", 1);

               # Numeric defaults to true because we are querying using dbID
               $list_of_features = $adaptor->_uncached_fetch_all_by_id_list(
                   [1,2,3,4],
                   undef,
                   "dbID");
  Description: This is a generic method used to fetch lists of features by IDs.
               It avoids caches, meaning it is best suited for block fetching.
               See fetch_all_by_dbID_list() for more info.
  Returntype : ArrayRef of Bio::EnsEMBL::Feature
  Exceptions : Thrown if a list of IDs is not supplied.
  Caller     : BaseFeatureAdaptor, BaseAdaptor and derived classes.

=cut

sub _uncached_fetch_all_by_id_list {
    my ( $self, $id_list_ref, $slice, $id_type, $numeric, $max_size ) = @_;

  if ( !defined($id_list_ref) || ref($id_list_ref) ne 'ARRAY' ) {
    throw("id_list list reference argument is required");
  }

  if ( !@{$id_list_ref} ) { return [] }

  # Construct a constraint like 't1.table1_id = 123'
  my @tabs = $self->_tables();
  my ( $name, $syn ) = @{ $tabs[0] };

  # prepare column name for query. If the id_type was
  # set to dbID then we assume the column must be
  # tablename_id e.g. gene_id. Otherwise we assume the id_type
  # is the field/column name
  my $field_name;
  if($id_type eq 'dbID') {
    $field_name = $name.'_id';
    # If numeric was not set default it to 1 since this is an int
    $numeric = 1 if ! defined $numeric;
  }
  elsif($id_type eq 'stable_id') {
    # If numeric was not set default it to 0 since this is a string
    $numeric = 0 if ! defined $numeric; 
    $field_name = $id_type;
  }
  else {
    $field_name = $id_type;
  }

  my $sql_data_type;

  # Ensuring we do not exceed MySQL's max_allowed_packet (defaults to 1MB)
  # by splitting large queries into smaller queries of at most 256KB
  # (262,144 8-bit characters)
  # If we had a numeric then really we are talking about working with
  # integers. Normal max ensembl id size is 12 plus 2 characters for
  # commas in our IN statement comes to 14. Even bloating this to 16 gives
  # a max number of 16,384 IDs (262114/16).
  # 
  if($numeric) {
    my $first_id = $id_list_ref->[0];
    if(!looks_like_number($first_id)) {
      throw "You specified that we are looking for numerics but $first_id is not a numeric";
    }
    $max_size = 16384 if ! defined $max_size;
    $sql_data_type = SQL_INTEGER;
  }
  # However when dealing with Strings those can be very large (assuming
  # 128 is the max length of a stable ID). 128 is 8x smaller than our
  # previous max expected integer so we reduce the max ids by 8. This gives
  # 2048 IDs (16384/8)
  else {
    $max_size = 2048 if ! defined $max_size;
    $sql_data_type = SQL_VARCHAR;
  }
  
  # build up unique id list, also validate on the way by
  my %id_list;
  for (@{$id_list_ref}) {
    $id_list{$_}++;
  }
  my @id_list = keys %id_list;

  my @out;
  my $inline = 1;
  while (@id_list) {
    my @ids;
    my $id_str;

    if ( scalar(@id_list) > $max_size ) {
      @ids = splice( @id_list, 0, $max_size );
    } 
    else {
      @ids     = @id_list;
      @id_list = ();
    }
    # Push off to our IN statement constructor for this work
    my $constraint = $self->generate_in_constraint(\@ids, "${syn}.${field_name}", $sql_data_type, $inline);
    push @out, @{ $self->generic_fetch($constraint, undef, $slice) };
  }

  return \@out;
}

# might not be a good idea, but for convenience
# shouldnt be called on the BIG tables though

sub fetch_all {
  my $self = shift;
  return $self->generic_fetch();
}

=head2 last_insert_id

  Arg [1]     : (optional) $field the name of the field the inserted ID was pushed 
                into
  Arg [2]     : (optional) HashRef used to pass extra attributes through to the 
                DBD driver
  Arg [3]     : (optional) $table the name of the table to use if the adaptor
                does not implement C<_tables()>
  Description : Delegating method which uses DBI to extract the last inserted 
                identifier. If using MySQL we just call the DBI method 
                L<DBI::last_insert_id()> since MySQL ignores any extra
                arguments. See L<DBI> for more information about this 
                delegated method. 
  Example     : my $id = $self->last_insert_id('my_id'); my $other_id = $self->last_insert_id();
  Returntype  : Scalar or undef
  
=cut

sub last_insert_id {
  my ($self, $field, $attributes, $table) = @_;
  my $dbc = $self->dbc();
  my $dbh = $dbc->db_handle();
  my @args;
  unless (@args = $dbc->_driver_object->last_insert_id_args($field, $table)) {
    if(!$table) {
      my ($table_entry) = $self->_tables(); # first table entry
      $table = $table_entry->[0];           # table_entry is [ name, alias ]
    }
    @args = (undef, $dbc->dbname(), $table, $field);
  }
  $attributes ||= {};
  return $dbh->last_insert_id(@args, $attributes);
}

=head2 insert_ignore_clause
=cut

sub insert_ignore_clause {
    my $self = shift;
    return $self->dbc->_driver_object->insert_ignore_clause;
}

=head2 _id_cache

  Description : Used to return an instance of a support BaseCache module
                which can be used to speed up object access. The method
                also respects the DBAdaptor's no_cache() flag and will
                return undef in those situations
  Example     : my $cache = $self->_id_cache();
  Returntype  : Bio::EnsEMBL::DBSQL::Support::BaseCache
  
=cut

sub _id_cache {
  my ($self) = @_;
  return if $self->db()->no_cache() && !$self->ignore_cache_override;
  if(! exists $self->{_id_cache}) {
    $self->{_id_cache} = $self->_build_id_cache();
  }
  return $self->{_id_cache};
}

=head2 _no_id_cache

  Description : Flags if the ID based caching is active or not. This could be
                due to the adaptor not wanting to cache or because of
                a global no_cache() flag on the DBAdaptor instance
  Returntype  : Boolean
  
=cut

sub _no_id_cache {
  my ($self) = @_;
  return 1 if ! $self->_id_cache();
  return 0;
}

=head2 ignore_cache_override

    Description : Method to interfere with no_cache directive from Registry on
                  a per adaptor basis. This method should be called after new()
                  in order to trigger the _build_id_cache at first query.                  
    Example     : $adaptor->ignore_cache_override(1);              
    Returntype  : Boolean

=cut

sub ignore_cache_override {
    my $self = shift;
    $self->{'_override'} = shift if(@_);
    unless (defined($self->{'_override'})) {return}
    return $self->{'_override'}; 
}

=head2 schema_version

    Description : Returns the schema version of the currently connected
                  DBAdaptor. The subroutine also caches this value so
                  repeated calls continue to be speedy.                  
    Example     : $adaptor->schema_version();            
    Returntype  : Integer

=cut

sub schema_version {
  my ($self) = @_;
  return $self->{_schema_version} if exists $self->{_schema_version};
  my $mc = $self->db()->get_MetaContainer();
  return $self->{_schema_version} = $mc->get_schema_version();
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

#_build_id_cache

#  Example    : my $id_cache = $self->_build_id_cache
#  Description: ABSTRACT PROTECTED
#               The subclass is responsible for returning an instance
#               of the Bio::EnsEMBL::DBSQL::Support::BaseCache
#               which can be used to speed up ID based fetch operations
#  Returntype : Instance of Bio::EnsEMBL::DBSQL::Support::BaseCache
#  Exceptions : Could be thrown by the implementing sub-class 
#  Caller     : BaseAdaptor::_id_cache
sub _build_id_cache {
  return;
}

sub dump_data {
  my $self = shift;
  my $data = shift;
  deprecate('This method is deprecated and will be removed in e91. Please use the get_all_attributes() and add_attributes() methods of DnaDnaAlignFeature instead. In the more general case, many feature types allow attributes to be stored as well');
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
    deprecate('This method is deprecated and will be removed in e91. Please use the get_all_attributes() and add_attributes() methods of DnaDnaAlignFeature instead. In the more general case, many feature types allow attributes to be stored as well');
    $data =~ s/\n|\r|\f|(\\\\)//g;
    return eval ($data); ## no critic
}

#
# Given a logic name and an existing constraint this will
# add an analysis table constraint to the feature.  Note that if no
# analysis_id exists in the columns of the primary table then no
# constraint is added at all
#
sub _logic_name_to_constraint {
  my $self = shift;
  my $constraint = shift;
  my $logic_name = shift;

  return $constraint if(!$logic_name);

  #make sure that an analysis_id exists in the primary table
  my ($prim_tab) = $self->_tables();
  my $prim_synonym = $prim_tab->[1];

  my $found_analysis=0;
  foreach my $col ($self->_columns) {
    my ($syn,$col_name) = split(/\./,$col);
    next if($syn ne $prim_synonym);
    if($col_name eq 'analysis_id') {
      $found_analysis = 1;
      last;
    }
  }

  if(!$found_analysis) {
    warning("This feature is not associated with an analysis.\n" .
            "Ignoring logic_name argument = [$logic_name].\n");
    return $constraint;
  }

  my $aa = $self->db->get_AnalysisAdaptor();
  my $an = $aa->fetch_by_logic_name($logic_name);

  if ( !defined($an) ) {
    return;
  }

  my $an_id = $an->dbID();

  $constraint .= ' AND' if($constraint);
  $constraint .= " ${prim_synonym}.analysis_id = $an_id";
  return $constraint;
}


1;

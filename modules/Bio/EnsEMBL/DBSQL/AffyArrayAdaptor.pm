#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::AffyArrayAdaptor
#
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AffyArrayAdaptor

=head1 SYNOPSIS

my $arrayAdaptor = $db->get_AffyArrayAdaptor();

my $array = $arrayAdaptor->fetch_by_name( 'some name' );
my $arrays = $arrayAdaptor->fetch_all();


=head1 DESCRIPTION

The AffyArrayAdaptor module stores and retrieves AffyArray objects in the database it is connected to.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::AffyArrayAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::AffyArray;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : my $array = $arrayAdaptor->fetch_by_name( 'U133_X3P' );
  Description: Retrieves an AffyArray object by name from the database. Its created
               from the entries in affy_array table
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub fetch_by_name {
    my $self = shift;
    my $name = shift;
    
    my $result = $self->generic_fetch( "aa.name = '$name'" );
    if( scalar @$result > 1 ) {
	warn( "AffyArray $name is not unique, only one result is returned" );
    } 

    return $result->[0];
}


=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::AffyArray $array
  Example    : none
  Description: This function is soley intended to lazy load attributes into empty 
               AffyArray objects. You should not need to call this.
  Returntype : none
  Exceptions : none
  Caller     : lazy load Array attributes into empty array from AffyArray object
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub fetch_attributes {
    my $self = shift;
    my $array = shift;

    my $tmp_array = $self->fetch_by_dbID( $array->dbID() );
    %$array = %$tmp_array;
}


=head2 _tablename

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub _tables {
  my $self = shift;
  
  return ['affy_array', 'aa'];
}


=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a list of columns to use for queries
  Returntype : list of strings
  Exceptions : none
  Caller     : internal
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub _columns {
  my $self = shift;

  return qw( aa.affy_array_id aa.parent_array_id aa.probe_setsize aa.name );

}


=head2 _objs_from_sth

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: PROTECTED implementation of superclass abstract method.
               creates SimpleFeatures from an executed DBI statement handle.
  Returntype : list reference to Bio::EnsEMBL::AffyFeature objects
  Exceptions : none
  Caller     : internal
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  my ( @result, $array_id, $parent_id, $setsize, $name );
 
  $sth->bind_columns( \$array_id, \$parent_id, \$setsize, \$name );

  while( $sth->fetch() ) {
      my $array = Bio::EnsEMBL::AffyArray->new
	  ( -dbID => $array_id,
	    -adaptor => $self,
	    -name => $name,
	    -setsize => $setsize
	  );
      push( @result, $array );
      if( $parent_id ) {
	  my $parent_array = Bio::EnsEMBL::AffyArray->new
	  ( -dbID => $parent_id,
	    -adaptor => $self,
	  );
	  $array->superset( $parent_array );
      }
  }
  return \@result;
}


=head2 store

  Arg [1..]  : list of Bio::EnsEMBL::AffyArray @arrays
  Example    : $arrayAdaptor->store( $array1, $array2, $array3 )
  Description: Stores te given AffyArray objects in the database. There is no check wether
               they already exist, so only call once per Array. Sets dbID and adaptor on the
               objects that it stores.
  Returntype : none
  Exceptions : none
  Caller     : affy feature calculating scripts
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut


sub store {
    my $self = shift;
    my @args = @_;
    
    for my $array ( @args ) {
	if( !$array->isa( "Bio::EnsEMBL::AffyArray" )) {
	    warn( "Can only store AffyArrays" );
	    next;
	}

	if( $array->dbID() && $array->adaptor() == $self ) {
	    # is already stored
	    next;
	}

	my $superset = $array->superset();
	if( defined $superset && ! $superset->dbID() ) {
	    $self->store( $superset );
	}

	my $sth = $self->prepare( "INSERT INTO affy_array".
				  "( name, probe_setsize, parent_array_id ) ".
				  "VALUES( ?,?,? )" );
	$sth->bind_param(1,$array->name(),SQL_VARCHAR);
	$sth->bind_param(2,$array->setsize(),SQL_INTEGER);
	if( defined $superset ) {
	    $sth->bind_param(3,$superset->dbID(), SQL_INTEGER);
	} else {
	    $sth->bind_param(3,undef);
	}
	$sth->execute();
	my $dbID = $sth->{'mysql_insertid'};
	$array->dbID( $dbID );
	$array->adaptor( $self );
    }
}

	    


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$simple_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("affy_array");
}

1;

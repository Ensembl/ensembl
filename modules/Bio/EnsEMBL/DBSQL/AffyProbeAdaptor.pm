#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::AffyProbeAdaptor
#
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AffyProbeAdaptor

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::AffyProbeAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::AffyProbe;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 store

  Arg [1]    : list of Bio::EnsEMBL::AffyFeatures @afs
               the affy features to store in the database
  Example    : $affy_feature_adaptor->store(@affy_feats);
  Description: Stores a list of affy feature objects in the database
  Returntype : none
  Exceptions : thrown if @afs is not defined, if any of the features do not
               have an attached slice or atached probe.
               or if any elements of @afs are not Bio::EnsEMBL::AffyFeatures 
  Caller     : general

=cut

sub store{
  my ($self,@probes) = @_;

  if( scalar(@probes) == 0 ) {
    throw("Must call store with list of AffyFeatures");
  }

  my $sth = $self->prepare
    ("INSERT INTO affy_probe ( affy_probe_id, affy_array_id, probeset, name )".
     "VALUES (?,?,?,?)" );

  my $db = $self->db();

 PROBE:
  foreach my $probe ( @probes ) {
      
      if( !ref $probe || !$probe->isa("Bio::EnsEMBL::AffyProbe") ) {
	  throw("Probe must be an Ensembl AffyProbe, " .
		"not a [".ref($probe)."]");
      }

      if($probe->is_stored($db)) {
	  warning("AffyProbe [".$probe->dbID."] is already stored" .
		  " in this database.");
	  next PROBE;
      }

      # all arrays in the probe need to be stored
      my $arrays = $probe->get_all_AffyArrays();
      my @stored_arrays;

      for my $array ( @$arrays ) {
	  if( defined $array->dbID() ) {
	      push( @stored_arrays, $array );
	  }
      }

      if( ! @stored_arrays ) {
	  warn( "Probes need attached arrays to be stored in database" );
	  next PROBE;
      }

      my $dbID;
      my $sth;

      for my $array ( @stored_arrays ) {
	  if( defined $dbID ) {
	      $sth = $self->prepare( "INSERT INTO affy_probe( affy_probe_id," .
				     "affy_array_id, name, probeset ) ".
				     "VALUES( ?,?,?,? )" );
	      $sth->execute( $dbID, $array->dbID(), $probe->get_probename( $array->name()),
			     $probe->probeset() );
	  } else {
	      $sth = $self->prepare( "INSERT INTO affy_probe( " .
				     "affy_array_id, name, probeset ) ".
				     "VALUES( ?,?,? )" );
	      $sth->execute( $array->dbID(), $probe->get_probename( $array->name()),
			     $probe->probeset() );
	      $dbID = $sth->{'mysql_insertid'};
	      $probe->dbID( $dbID );
	      $probe->adaptor($self);
	  }
      }
  }
}

sub fetch_by_probeset {
    my $self = shift;
    my $probeset = shift;
    
    return $self->generic_fetch( "ap.probeset = $probeset" );
}



sub fetch_by_AffyArray {
    my $self = shift;
    my $array = shift;

    if( !ref( $array) || !$array->isa( "Bio::EnsEMBL::AffyArray" ) ) {
	warn( "Argument must be a stored AffyArray" );
	return [];
    }
    my $array_id = $array->dbID();
    if( ! defined $array_id ) {
	warn( "Argument must be a stored AffyArray" );
	return [];	
    }

    return $self->generic_fetch( "ap.affy_array_id = array_id" );
}

sub fetch_by_AffyFeature {
    my $self = shift;
    my $feature = shift;

    if( ! ref( $feature ) || !$feature->isa( "Bio::EnsEMBL::AffyFeature" ) ||
	! $feature->{'_probe_id'} ) {
	throw( "Need AffyFeature from database as argument" );
    }

    $self->fetch_by_dbID( $feature->{'_probe_id'} );
}



=head2 _tablename

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal

=cut

sub _tables {
  my $self = shift;
  
  return ['affy_probe', 'ap'];
}


=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a list of columns to use for queries
  Returntype : list of strings
  Exceptions : none
  Caller     : internal

=cut

sub _columns {
  my $self = shift;

  return qw( ap.affy_probe_id ap.affy_array_id ap.probeset ap.name );

}


=head2 _objs_from_sth

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: PROTECTED implementation of superclass abstract method.
               creates SimpleFeatures from an executed DBI statement handle.
  Returntype : list reference to Bio::EnsEMBL::AffyFeature objects
  Exceptions : none
  Caller     : internal

=cut

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  my ( @result, $current_dbid, $probe_id, $array_id, $probeset, $name );
  my ( $array, %array_cache) ;

  $sth->bind_columns( \$probe_id, \$array_id, \$probeset, \$name );
  my $probe;

  while( $sth->fetch() ) {
      $array  = ( $array_cache{$array_id} ||= 
		  $self->db->get_AffyArrayAdaptor()->fetch_by_dbID( $array_id ) );
      if( $current_dbid != $probe_id ) {
	  # make a new probe
	  $probe = Bio::EnsEMBL::AffyProbe->new(
	      -array => $array,
	      -probeset => $probeset,
	      -name => $name,
	  );
	  push( @result, $probe );
	  $current_dbid = $probe_id;
      } else {
	  $probe->add_Array_probename( $array, $name );
      }
  }
  return \@result;
}
 


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$simple_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("affy_probe");
}

1;

#
# Ensembl module for Bio::EnsEMBL::AffyArray
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AffyArray - A module to represent the name of an Affymetrix oligo array with
some attached functionality.

=head1 SYNOPSIS

use Bio::EnsEMBL::AffyArray;

$feature = Bio::EnsEMBL::AffyArray->new 
    ( -name => 'Affy-1',
      -setsize => 11,
      -included_in => $another_array.
      -probecount => 40000
    )      

=head1 DESCRIPTION

AffyArray contains the informmation in the affy_array table, currently the name, 
the set_size of the array, and if there is superset array.

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;


package Bio::EnsEMBL::AffyArray;

use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use vars qw(@ISA);

use Bio::EnsEMBL::Storable;

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg  [NAME] : 
               Name of the array
  Arg  [INCLUDED_IN]  : 
               A possible superset Array
  Arg  [SETSIZE]  : 
               How many probes per normal set on this Array
  Arg  [PROBECOUNT]  : 
               How many probes are on this Array
  Example    : none
  Description: 
  Returntype : 
  Exceptions : none
  Caller     : 
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub new {
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my($name, $superset, $probecount, $setsize ) =
    rearrange(['NAME', 'INCLUDED_IN', 'PROBECOUNT', 'SETSIZE'],
              @_);

  $self->{'name'} = $name if(defined $name);
  $self->{'superset'} = $superset if( defined $superset);

  $self->{'probecount'}    = $probecount if( defined  $probecount);
  $self->{'setsize'} = $setsize if( defined $setsize );

  return $self;
}


=head2 get_all_AffyProbes

  Args       : none 
  Example    : my $probes = $array->get_all_AffyProbes();
  Description: Returns all probes that are associated with this Array.
               Only works when connected with the database.
  Returntype : listref of Bio::EnsEMBL::AffyProbe
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub get_all_AffyProbes {
  my $self = shift;

  if( $self->adaptor() && $self->dbID() ) {
    my $probeAdaptor = $self->adaptor()->db()->get_AffyProbeAdaptor();
    my $probes = $probeAdaptor->fetch_by_AffyArray( $self );
    return $probes;
  } else {
    warning( "Need database connection to retrieve Probes" );
    return [];
  }
}

    

=head2 name

  Arg [1]    : string $name
  Example    : none
  Description: getter / setter / lazy load of attribute name of the
               AffyArray. The only mandatory attribute for ths object.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub name {
    my $self = shift;
    $self->{'name'} = shift if( @_ );
    if(( ! exists $self->{'name'}) && $self->{'dbID'} && $self->{'adaptor'} ) {
	$self->adaptor->fetch_attributes( $self );
    }
    return $self->{'name'};
}


=head2 setsize

  Arg [1]    : int $setsize
  Example    : none
  Description: getter / setter / lazy load of attribute setsize.
               The setsize is the number of probes making up one standard probeset
               on this Array. 
  Returntype : listref of Bio::EnsEMBL::AffyFeature
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub setsize {
    my $self = shift;
    $self->{'setsize'} = shift if( @_ );
    if(( ! exists $self->{'setsize'}) && $self->{'dbID'} && $self->{'adaptor'} ) {
	$self->adaptor->fetch_attributes( $self );
    }
    return $self->{'setsize'};
}




=head2 superset

  Arg [1]    : Bio::EnsEMBL::AffyArray $superset
  Example    : none
  Description: getter / setter / lazy load for attribute superset. A superset
               is another AffyArray that contains all the probes of this Array.
               This is entirely optional / (bordering to superfluous). 
  Returntype : Bio::EnsEMBL::AffyArray
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub superset {
    my $self = shift;
    if( @_ ) {
	my $superset = shift;
	if( $superset ) {
	    throw( "Need Bio::EnsEMBL::AffyArray as superset" ) unless
		(ref( $superset ) && $superset->isa( "Bio::EnsEMBL::AffyArray" ));
	}
	$self->{'superset'} = $superset;
    } elsif (( ! exists $self->{'superset'}) && $self->{'dbID'} && $self->{'adaptor'} ) {
	$self->adaptor->fetch_attributes( $self );
    }
    return $self->{'superset'};
}








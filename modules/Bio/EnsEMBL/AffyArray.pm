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

AffyArray contains the informmation in the affy_array table, currently the name, the set_size of the array,
and if there is superset array.

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;


package Bio::EnsEMBL::AffyArray;

use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Utils::Exception qw( throw );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);

sub new {
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my($name, $superset, $probecount, $setsize ) =
    rearrange(['NAME', 'INCLUDED_IN', 'PROBECOUNT', 'SETSIZE'],
              @_);

  $self->{'name'} = $name;
  $self->{'superset'} = $superset;

  $self->{'probecount'}    = $probecount;
  $self->{'setsize'} = $setsize;

  return $self;
}


    
sub name {
    my $self = shift;
    $self->{'name'} = shift if( @_ );
    return $self->{'name'};
}

sub setsize {
    my $self = shift;
    $self->{'setsize'} = shift if( @_ );
    return $self->{'setsize'};
}

sub probecount {
    my $self = shift;
    $self->{'probecount'} = shift if( @_) ;
    return $self->{'probecount'};
}

sub superset {
    my $self = shift;
    if( @_ ) {
	my $superset = shift;
	if( $superset ) {
	    throw( "Need Bio::EnsEMBL::AffyArray as superset" ) unless
		(ref( $superset ) && $superset->isa( "Bio::EnsEMBL::AffyArray" ));
	}
	$self->{'superset'} = $superset;
    }
    return $self->{'superset'};
}








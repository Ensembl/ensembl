#
# Ensembl module for Bio::EnsEMBL::AffyProbe
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AffyProbe - a module to represent the data for an affymetrix oligo spot.
The same oligo as part of the same Set can be found on many Arrays.


=head1 SYNOPSIS

use Bio::EnsEMBL::AffyProbe;

$feature = Bio::EnsEMBL::AffyProbe->new 
    
      -probeset => 'some setname'
)

=head1 DESCRIPTION

Affyprobe represent an oligo (probe) on one or more Affymetrix Arrays. 

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::AffyProbe;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);

sub new {
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my($arrays, $probenames, $probeset, $arraynames ) =
      rearrange(['ARRAYS', 'PROBENAMES', 'PROBESET', 'ARRAYNAMES' ],
		@_);

  $self->arrays( $arrays ) if defined( $arrays );
  $self->probenames( $probenames ) if defined $probenames;
  $self->arraynames( $arraynames ) if $arraynames;


  $self->{'probeset'} = $probeset;

  return $self;
}


sub add_Array {
    my $self = shift;
    my ( $array, $probename ) = @_;
}

sub get_all_AffyFeatures {
    my $self = shift;
}

sub get_all_AffyArrays {
    my $self = shift;
}

sub get_all_complete_names {
    my $self = shift;
}

sub get_complete_name {
    my $self = shift;
    my $arrayname = shift;
}

sub get_all_probenames {
    my $self = shift;
}

sub get_probename {
    my $self = shift;
    my $arrayname = shift;
}

sub add_array_name {
    my $self = shift;
    my ( $affy_array, $probename ) = @_;
}



# if the probe is set, take it from there
sub probeset {
    my $self = shift;
    $self->{'probeset'} = shift if( @_ );
    return $self->{'probeset'};
}

# lazy load this
sub probe {
    my $self = shift;
    if( @_ ) {
	my $probe = shift;
	if( $probe ) {
	    throw( "Need Bio::EnsEMBL::AffyProbe as probe" ) unless
		(ref( $probe ) && $probe->isa( "Bio::EnsEMBL::AffyProbe" ));
	}
	$self->{'probe'} = $probe;
    } else {
	if( defined $self->adaptor() && $self->dbID() ) {
	    $self->{'probe'} = $self->adaptor()->db()->
		get_AffyProbeAdaptor()->fetch_by_AffyFeature( $self );
	}
    }
    return $self->{'probe'};
}








#
# Ensembl module for Bio::EnsEMBL::DensityFeature
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DensityFeature - A feature representing a density, or precentage
coverage etc. in a given region.

=head1 SYNOPSIS

use Bio::EnsEMBL::DensityFeature;

$feature = Bio::EnsEMBL::DensityFeature->new(-start    => 100,
                                             -end      => 220,
                                             -slice    => $slice,
                                             -analysis => $analysis,
                                             -density_value    => 90.1,
                                             -dbID     => 112,
                                             -adaptor  => $adaptor);

=head1 DESCRIPTION

A density feature represents a count, density, or percentage coverage, etc. for
a given region.

This module is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::DensityFeature;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [DENSITY_VALUE] : The number of features which were found within the
               region of this DensityFeature.  May also be a percentage or
               coverage, etc.
  Arg [...]  : Named arguments passed to superclass
  Example    : $feature = Bio::EnsEMBL::SimpleFeature->new
                            (-start    => 1,
                             -end      => 1e6,
                             -analysis => $analysis,
                             -density_value => 80.5);
  Description: Creates a new density feature.
  Returntype : Bio::EnsEMBL::DensityFeature
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;

 #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my($density_value) = rearrange(['DENSITY_VALUE'], @_);

  throw("Density value must be >= 0.") if($density_value < 0);

  $self->{'density_value'} = $density_value;
  $self->{'strand'} = 0;

  return $self;
}


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless($hashref,$class);
}


=head2 strand

  Arg [1]    : none
  Example    : $strand = $df->strand();
  Description: Getter fot the strand attribute. Density features always have
               strand 0 and this attribute is not settable.
  Returntype : int (always 0)
  Exceptions : warning if an attempt is made to set the strand
  Caller     : general

=cut

sub strand {
  my $self = shift;
  warning("DensityFeature strand is not settable") if(@_);
  return 0;
}



=head2 density_value

  Arg [1]    : (optional) float $density_value
  Example    : $dv = $density_feature->density_value();
  Description: Getter/Setter for the density value of this DensityFeature.
               The density value may be a count, a percentage, or a coverage
               of a feature type in the area defined by this feature.
  Returntype : float
  Exceptions : throw if a negative density value is provided
  Caller     : general

=cut

sub density_value {
  my $self = shift;

  if(@_) {
    my $density_value = shift;
    throw("Density value must be >= 0.") if($density_value < 0);
    $self->{'density_value'} = $density_value;
  }

  return $self->{'density_value'};
}


1;




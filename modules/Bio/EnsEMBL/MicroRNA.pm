=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::MicroRNA - A class representing a microRNA product
of a transcript

=head1 DESCRIPTION

A specialisation of Bio::EnsEMBL::RNAProduct describing
MicroRNAs. Mostly takes care of wrapping miRNA-specific RNAProduct
attributes in methods which make them look like ordinary class
members.

=head1 SYNOPSIS

  my $miR = Bio::EnsEMBL::MicroRNA->new(
    -SEQ_START => 36,
    -SEQ_END   => 58
  );

  # Stable-ID setter
  $miR->stable_id('ENSS00090210');

  # Get start and end position in the precursor transcript
  my $start = $miR->start();
  my $end = $miR->end();

=cut


package Bio::EnsEMBL::MicroRNA;

use vars qw($AUTOLOAD);
use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref wrap_array );
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::RNAProduct;

use parent qw(Bio::EnsEMBL::RNAProduct);


=head2 new

  Arg: [-ARM]         : which arm of the hairpin precursor this miRNA comes
                        from. Returns 3 and 5 for 3' and 5', respectively.
  Arg [...]           : Named arguments to superclass constructor
                        (see Bio::EnsEMBL::RNAProduct)
  Example    : my $miR = Bio::EnsEMBL::MicroRNA->new(
                 -SEQ_START => 36,
                 -SEQ_END   => 58,
                 -ARM       => 3
               );
  Description: Constructor.  Creates a new MicroRNA object
  Returntype : Bio::EnsEMBL::MicroRNA
  Exceptions : throw if ARM value is out of bounds
  Caller     : general
  Status     : In Development

=cut

# perlcritic doesn't know about rearrange(), silence it
sub new { ## no critic (Subroutines::RequireArgUnpacking)
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($arm) = rearrange(["ARM"], @_);
  if (defined($arm)) {
    _validate_arm_value($arm);
  }
  $self->{'arm'} = $arm;

  return $self;
}


=head2 arm

    Arg [1]     : (optional) int $arm which arm of the hairpin precursor
                  this miRNA comes from
    Example     : $mirna_arm = $mirna->arm();
                  $mirna->arm(3);
    Description : Sets or returns the arm of the hairpin this miRNA comes
                  from. Accepted values are 3 and 5 for 3' and 5',
                  respectively.
    Return type : Integer
    Exceptions  : throw if setter is passed an incorrect value
                  or if multiple 'mirna_arm' attributes exist.
    Caller      : General
    Status      : Stable

=cut

sub arm {
  my ($self, $arm) = @_;

  if (defined $arm) {
    _validate_arm_value($arm);
    $self->{'arm'} = $arm;
  } elsif (!defined($self->{'arm'})) {
    my $arm_attrs = $self->get_all_Attributes('mirna_arm');
    my $n_arms = scalar @{$arm_attrs};
    if ($n_arms > 0) {
      if ($n_arms > 1) {
        throw("MicroRNA " . $self->display_id() .
              " has multiple arm attributes");
      }
      $self->{'arm'} = $arm_attrs->[0]->value();
    }
  }

  return $self->{'arm'};
}


=head2 summary_as_hash

  Example       : $mirna_summary = $mirna->summary_as_hash();
  Description   : Retrieves a textual summary of this MicroRNA.
                  Built on top of generic implementation in RNAProduct.
  Returns       : hashref of arrays of descriptive strings
  Status        : Intended for internal use

=cut

sub summary_as_hash {
  my $self = shift;

  my $summary = SUPER::summary_as_hash();
  $summary->{'arm'} = $self->arm();

  return $summary;
}


=head2 _validate_arm_value
  Arg [1]    : int $arm which arm of the hairpin precursor this miRNA
               comes from
  Description: PRIVATE validates if its argument has one of the accepted
               values for specifying the miRNA hairpin arm.
  Returntype : none
  Exceptions : throw if the argument is out of bounds
  Caller     : internal
  Status     : Stable

=cut

sub _validate_arm_value {
  my ($arm) = @_;

  if (($arm != 3) && ($arm != 5)) {
    throw("'$arm' is not a valid miRNA hairpin-arm specification");
  }

  return;
}

1;

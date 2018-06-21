=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

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

TODO

=head1 SYNOPSIS

  my $miR = Bio::EnsEMBL::MicroRNA->new(
    -SEQ_START => 36,
    -SEQ_END   => 58
  );

  # Stable-ID setter
  $miR->stable_id('ENSM00090210');

  # Get start and end position in the precursor transcript
  my $start = $miR->start();
  my $end = $miR->end();

=cut


package Bio::EnsEMBL::MicroRNA;

use vars qw($AUTOLOAD);
use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref wrap_array );
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::RNAProduct;

use parent qw(Bio::EnsEMBL::RNAProduct);


=head2 new

  Arg [-SEQ_START]    : The offset in the Translation indicating the start
                        position of the product sequence.
  Arg [-SEQ_END]      : The offset in the Translation indicating the end
                        position of the product sequence.
  Arg [-STABLE_ID]    : The stable identifier for this RNAPRoduct
  Arg [-VERSION]      : The version of the stable identifier
  Arg [-DBID]         : The internal identifier of this MicroRNA
  Arg [-ADAPTOR]      : The TranslationAdaptor for this MicroRNA
  Arg [-SEQ]          : Manually sets the nucleotide sequence of this
                        rnaproduct. May be useful if this rnaproduct is not
                        stored in a database.
  Arg [-CREATED_DATE] : the date the rnaproduct was created
  Arg [-MODIFIED_DATE]: the date the rnaproduct was modified
  # FIXME: does all of the above have to be repeated here or can we somehow pull this from the superclass?
  Arg: [-ARM]         : which arm of the hairpin precursor this miRNA comes
                        from. Negative values indicate 3', positive ones - 5'.
                        FIXME: is this right?
  Example    : my $miR = Bio::EnsEMBL::MicroRNA->new(
                 -SEQ_START => 36,
                 -SEQ_END   => 58,
                 -ARM       => -1
               );
  Description: Constructor.  Creates a new MicroRNA object
  Returntype : Bio::EnsEMBL::MicroRNA
  Exceptions : none
  Caller     : general
  Status     : In Development

=cut

sub new { ## no critic (Subroutines::RequireArgUnpacking)
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  # FIXME: see the comment about same in RNAProduct::new()
  my $type_id = 1;

  my ($arm) = rearrange(["ARM"], @_);
  $self->{'arm'} = $arm;

  return $self;
}


=head2 arm

    Example     : $mirna_arm = $mirna->arm();
    Description : Returns the arm of the hairpin this miRNA comes from.
                  Negative values indicate the 3' end, positive ones
                  - the 5' one. FIXME: is this right?
    Return type : Integer
    Exceptions  : warn if multiple 'mirna_arm' attributes exist
    Caller      : General
    Status      : Stable

=cut

sub arm {
  my ($self, $arm) = @_;

  if (defined $arm) {
    $self->{'arm'} = $arm;
  } elsif (!defined($self->{'arm'})) {
    my $arm_attrs = $self->get_all_Attributes('mirna_arm');
    my $n_arms = scalar @{$arm_attrs};
    if ($n_arms > 0) {
      if ($n_arms > 1) {
	carp("MicroRNA " . $self->display_id() .
	     " has multiple arm attributes, using first");
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

1;

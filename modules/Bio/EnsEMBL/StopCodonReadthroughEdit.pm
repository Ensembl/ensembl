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

  Bio::EnsEMBL::StopCodonReadthroughEdit - Object representing a stop codon readthrough edit in a sequence

=head1 SYNOPSIS

  use Bio::EnsEMBL::StopCodonReadthroughEdit;

  # Get transcript
  my $transcript_adaptor = $db->get_TranscriptAdaptor();
  my $transcript = $transcript_adaptor->fetch_by_stable_id("ENST00000217347");
  $transcript->edits_enabled(1);
  print "Before modifiction: $transcript->translate->seq()\n";

  # Construct a stop codon readthrough edit object
  my $stop_codon_readthrough_edit = Bio::EnsEMBL::StopCodonReadthroughEdit->new(265);

  # Apply post translation edit
  $transcript->translation->add_Attributes($stop_codon_readthrough_edit->get_Attribute());
  my $translated_sequence = $transcript->translate->seq();
  print "After modifiction: $translated_sequence\n";

=head1 DESCRIPTION

  Biologically, STOP codon readthrough is a rare phenomenon whereby translation
  does not terminate at an in-frame STOP codon, but instead continues further downstream.
  It is believed the STOP codon is instead read as a 'sense' codon, i.e. encodes for an amino acid.

  The location of a STOP codon readthrough is indicated by an asterisk (*) in the post translation sequence.
  This class edits the sequence to replace the asterisk (*) with an 'X' to make it similar to RefSeq representation.

=head1 METHODS

=cut

package Bio::EnsEMBL::StopCodonReadthroughEdit;

use strict;
use warnings;

use parent qw(Bio::EnsEMBL::SeqEdit);

=head2 new

  Arg [-POSITION]  :
       int - start and end position of the stop codon readthrough edit in the sequence

  Example    : $stop_codon_rt_edit = Bio::EnsEMBL::StopCodonReadthroughEdit->new($position);
  Description: Creates a new stop codon readthrough edit object
  Returntype : Bio::EnsEMBL::StopCodonReadthroughEdit
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ($self, $position) = @_;

  my $class = ref($self) || $self;
  my $stop_codon_rt_edit = $class->SUPER::new(
        -START   => $position,
        -END     => $position,
        -ALT_SEQ => 'X',
        -CODE    => '_stop_codon_rt');

  return $stop_codon_rt_edit;

}

1;

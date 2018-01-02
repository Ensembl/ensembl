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

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::VegaCuration::Gene;

use strict;
use warnings;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::ConversionSupport;

@ISA = qw(Bio::EnsEMBL::Utils::ConversionSupport);


=head2 find_gaps

   Args       : arrayref of B::E::Transcripts
   Example    : my $gaps = find_gaps($all_transcripts)
   Description: identifies regions of a gene that are not covered by any transcript
   Returntype : int
   Exceptions : none
   Caller     : internal

=cut

sub find_gaps {
  my $self = shift;
  my ($all_transcripts) = @_;
  my $gaps = 0;
  my @sorted_transcripts = sort {$a->start <=> $b->start || $b->end <=> $a->end} @{$all_transcripts};
  if ( my $first_transcript = shift @sorted_transcripts ) {
    my $pos = $first_transcript->end;
    foreach my $transcript (@sorted_transcripts) {
      next if ($transcript->end < $pos );
      if ($transcript->start < $pos && $transcript->end > $pos ) {
	$pos = $transcript->end;			
	next;
      }
      elsif ($transcript->end > $pos) {
	$gaps++;
	$pos = $transcript->end;
      }
    }
  }
  return $gaps;		
}

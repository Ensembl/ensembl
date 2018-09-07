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

package XrefParser::WormbaseCElegansUniProtParser;

# UniProt xrefs are sometimes - really - dependent xrefs of
# INSDC entries which we get from somewhere else
# Attempt to find the parent (has to already be present in the xref table)
# INSDC and UniProt entries have the same protein sequence, and
# UniProt lists INSDC as a parent. We get INSDC entries from somewhere else,
# so make UniProt entries dependent on INSDC entries.
# Note:
# INSDC entries have coordinates, and UniProt entries don't.
# So for perfect homologs, there can be many INSDC entries per UniProt.

use parent XrefParser::WormbaseCElegansBase, XrefParser::UniProtParser;

sub upload_xref_object_graphs {
  my ($self, $xrefs, $dbi) = @_;
  my $source_id = $self->get_source_id_for_source_name('protein_id'); 
  my $source_id_skip = $self->get_source_id_for_source_name('EMBL'); 
  my @adapted_xrefs;
  for my $xref ( @$xrefs) {
    push @adapted_xrefs, $self->swap_dependency($source_id, $dbi, $xref, $source_id_skip);
  }  
  return $self->SUPER::upload_xref_object_graphs(\@adapted_xrefs, $dbi);
}
1;

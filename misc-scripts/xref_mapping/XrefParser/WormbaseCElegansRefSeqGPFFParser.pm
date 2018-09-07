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

package XrefParser::WormbaseCElegansRefSeqGPFFParser;

use parent XrefParser::WormbaseCElegansBase, XrefParser::RefSeqGPFFParser;

my $source_id;

sub upload_xref_object_graphs {
  my ($self, $xrefs, $dbi) = @_;
  $source_id //= $self->get_source_id_for_source_name('protein_id'); 
  my @adapted_xrefs;
  for my $xref ( @$xrefs) {
    push @adapted_xrefs, $self->swap_dependency($source_id, $dbi, $xref);
  }  
  return $self->SUPER::upload_xref_object_graphs(\@adapted_xrefs, $dbi);
}
sub xref_from_record {
   my ($self, $entry, @args) = @_;
   
   my $xref = $self->SUPER::xref_from_record($entry, @args);
   $source_id //= $self->get_source_id_for_source_name('protein_id'); 
   $entry =~ /This record has been curated by WormBase. The\s+reference sequence is identical to (.*?)\./;
   my $insdc_protein_id = $1;
   if($insdc_protein_id) {
     $xref->{DEPENDENT_XREFS} //= [];
     push @{$xref->{DEPENDENT_XREFS}}, {ACCESSION => $insdc_protein_id, SOURCE_ID=>$source_id};
     return $xref;
   } else {
     return undef;
   }
}
1;

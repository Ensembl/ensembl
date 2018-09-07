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

package XrefParser::WormbaseCElegansBase;

sub swap_dependency {
  my ($self, $source_id, $dbi, $xref, @source_ids_skip) = @_;

  my @matching_source_id_dependents;
  my @other_dependents;
  for my $dependent_xref (@{$xref->{DEPENDENT_XREFS} || []}){
     my $source_id_here = $dependent_xref->{SOURCE_ID};
     if($source_id_here eq $source_id
         and $self->get_xref($dependent_xref->{ACCESSION}, $dependent_xref->{SOURCE_ID}, $xref->{SPECIES_ID})){
         $dependent_xref->{SPECIES_ID} = $xref->{SPECIES_ID};
         push @matching_source_id_dependents, $dependent_xref;
     } elsif (grep {$_ == $source_id_here} @source_ids_skip){
       #skip
     } else {
         push @other_dependents, $dependent_xref;
     }
  }
  return map {{%$_, LABEL=>undef, INFO_TYPE => "MISC", DEPENDENT_XREFS => [{
        %$xref,
        INFO_TYPE => "DEPENDENT",
        LINKAGE_SOURCE_ID => $source_id,
     }, map {{%$_,INFO_TYPE => "DEPENDENT", LINKAGE_SOURCE_ID => $source_id}} @other_dependents]}} @matching_source_id_dependents;
}
1;

=head1 LICENSE

Copyright [2018-2020] EMBL-European Bioinformatics Institute

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
use strict;
use warnings;
package XrefParser::WormbaseCElegansBase;

sub swap_dependency {
  my ($self, $source_ids, $dbi, $xref, @source_ids_skip) = @_;
  $source_ids = [$source_ids] unless ref $source_ids eq 'ARRAY';

  my @matching_source_id_dependents;
  my @other_dependents;
  for my $dependent_xref (@{$xref->{DEPENDENT_XREFS} || []}){
     my $source_id_here = $dependent_xref->{SOURCE_ID};
     if(grep {$_ == $source_id_here } @$source_ids
         and $self->get_xref($dependent_xref->{ACCESSION}, $dependent_xref->{SOURCE_ID}, $xref->{SPECIES_ID})){
         $dependent_xref->{SPECIES_ID} = $xref->{SPECIES_ID};
         push @matching_source_id_dependents, $dependent_xref;
     } elsif (grep {$_ == $source_id_here} @source_ids_skip){
       #skip
     } else {
         push @other_dependents, $dependent_xref;
     }
  }
  my @result;
  for my $matching_source_id_dependent (@matching_source_id_dependents) {
    my $source_id = $matching_source_id_dependent->{SOURCE_ID};
    my $xref_as_dependent_here = {
      %$xref, INFO_TYPE => "DEPENDENT",
      LINKAGE_SOURCE_ID => $source_id,
      DEPENDENT_XREFS => undef,
    };
    my @dependents_here = ({
      %$xref, INFO_TYPE => "DEPENDENT",
      LINKAGE_SOURCE_ID => $source_id,
      DEPENDENT_XREFS => undef,
    });
    for my $d (@other_dependents){
       push @dependents_here, {
          %$d, INFO_TYPE => "DEPENDENT", LINKAGE_SOURCE_ID => $source_id,
       };
    }
    push @result, {
       %$matching_source_id_dependent,
       LABEL=>undef, INFO_TYPE => "MISC",
       DEPENDENT_XREFS =>\@dependents_here, 
    };
  }

  return @result;
}
1;

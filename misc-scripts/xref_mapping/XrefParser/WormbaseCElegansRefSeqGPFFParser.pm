=head1 LICENSE

Copyright [2018-2019] EMBL-European Bioinformatics Institute

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
use strict;
use warnings;
use parent qw/XrefParser::WormbaseCElegansBase XrefParser::RefSeqGPFFParser/;
my $SOURCE_IDS;
my $ACCESSION_FROM_ENTRY_PATTERN;
sub run {
  my ($self, $arg_ref) = @_;
  my $type = $self->type_from_file(@{$arg_ref->{files}});
  if($type eq 'peptide'){
    $SOURCE_IDS = [ $self->get_source_id_for_source_name('protein_id') ];
    $ACCESSION_FROM_ENTRY_PATTERN =  qr/This record has been curated by WormBase. The\s+reference sequence is identical to (.*?)\./;
  } elsif ($type eq 'dna'){
    $SOURCE_IDS = [
       $self->get_source_id_for_source_name('wormbase_cds'),
       $self->get_source_id_for_source_name('wormbase_transcript'),
    ];
    $ACCESSION_FROM_ENTRY_PATTERN = qr/standard_name="(.*?)"/;
  }
  die %$arg_ref unless @$SOURCE_IDS;
  return $self->SUPER::run($arg_ref);
}
sub upload_xref_object_graphs {
  my ($self, $xrefs, $dbi) = @_;
  my @adapted_xrefs;
  for my $xref ( @$xrefs) {
    push @adapted_xrefs, $self->swap_dependency($SOURCE_IDS, $dbi, $xref, @$SOURCE_IDS);
  }  
  return $self->SUPER::upload_xref_object_graphs(\@adapted_xrefs, $dbi);
}

sub xref_from_record {
   my ($self, $entry, @args) = @_;
   return &modify_xref_with_dependent(
      $SOURCE_IDS, $entry,
      $self->SUPER::xref_from_record($entry, @args),
      $ACCESSION_FROM_ENTRY_PATTERN,
   );
}

sub modify_xref_with_dependent {
   my ($source_ids, $entry, $xref, $get_accession_pattern) = @_;
   return unless $xref;
   my ($accession) = $entry =~ $get_accession_pattern;
   return unless $accession;
   $xref->{DEPENDENT_XREFS} //= [];
   push @{$xref->{DEPENDENT_XREFS}}, map {{ACCESSION => $accession, SOURCE_ID=>$_}} @$source_ids;
   return $xref;
}
1;

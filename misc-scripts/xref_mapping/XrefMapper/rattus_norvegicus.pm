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

package XrefMapper::rattus_norvegicus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_official_name {
  return 'RGD';
}

sub get_canonical_name{
   return "RGD";
}

# Not running transcript_names_from_gene for merged species
# as this is already been done in the OfficialNaming mapper
sub transcript_names_from_gene {
  return;
}

sub set_display_xrefs{
  my $self = shift;
  my $display = XrefMapper::DisplayXrefs->new($self);
  $display->set_display_xrefs_from_stable_table();

}

sub set_gene_descriptions(){
  my $self = shift;
  my $display = XrefMapper::DisplayXrefs->new($self);
  $display->set_gene_descriptions_from_display_xref()
}

sub gene_description_filter_regexps {

  return ('^\(CLONE REM\d+\) ORF \(FRAGMENT\)\.*',
	  '^ORF\s*\d+\s+PROTEIN\.*',
	  '\(?[0-9A-Z]{10}RIK PROTEIN\)?[ \.]',
	  'RIKEN CDNA [0-9A-Z]{10}[ \.;]',
	  '.*RIKEN FULL-LENGTH ENRICHED LIBRARY.*PRODUCT:',
	  '.*RIKEN FULL-LENGTH ENRICHED LIBRARY.*',
	  '\(*HYPOTHETICAL\s+.*',
	  '^UNKNOWN\s+.*',
	  'CDNA SEQUENCE\s?,? [A-Z]+\d+[ \.;]',
	  'CLONE MGC:\d+[ \.;]',
	  'MGC:\s*\d+[ \.;]',
	  'HYPOTHETICAL PROTEIN,',
	  'HYPOTHETICAL PROTEIN \S+[\.;]',
	  'DNA SEGMENT, CHR.*',
	  'PROTEIN \S+ HOMOLOG\.?',
	  '^SIMILAR TO GENE.*',
	  'SIMILAR TO PUTATIVE[ \.]',
	  '^SIMILAR TO HYPOTHETICAL.*',
	  'SIMILAR TO (KIAA|LOC|RIKEN).*',
	  'SIMILAR TO GENBANK ACCESSION NUMBER\s+\S+',
	  'SIMILAR TO\s+$',
          'EXPRESSED SEQUENCE [A-Z]+\d+[ \.;]',
          'EST [A-Z]+\d+[ \.;]',
          '^\s*\(FRAGMENT\)\.?\s*$',
	  '^\s*\(?GENE\)?\.?;?\s*$',
          '\s*\(?GENE\)?\.?;?',
          '\s*\(?PRECURSOR\)?\.?;?',
          '^\s*\(\s*\)\s*$',
	  '^\s*\(\d*\)\s*[ \.]$',
          '^\s+\(?\s*$');

}

1;

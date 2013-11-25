=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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


sub gene_display_xref_sources {
  my $self     = shift;
	
  my @list = qw(RFAM
                miRBase
                RGD
                Uniprot_gn
                EntrezGene);

  my %ignore;

  #don't use EntrezGene labels dependent on predicted RefSeqs

$ignore{'EntrezGene'} =<<IEG;
SELECT DISTINCT ox.object_xref_id
  FROM object_xref ox, dependent_xref dx, 
       xref xmas, xref xdep, 
       source smas, source sdep
    WHERE ox.xref_id = dx.dependent_xref_id AND
          dx.dependent_xref_id = xdep.xref_id AND
          dx.master_xref_id = xmas.xref_id AND
          xmas.source_id = smas.source_id AND
          xdep.source_id = sdep.source_id AND
          smas.name like "Refseq%predicted" AND
          sdep.name like "EntrezGene" AND
          ox.ox_status = "DUMP_OUT" 	 
IEG

  #don't use labels starting with LOC

$ignore{'LOC_prefix'} =<<LOCP;
SELECT object_xref_id
  FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
   WHERE ox_status = 'DUMP_OUT' AND label REGEXP '^LOC[[:digit:]]+'
LOCP

  return [\@list,\%ignore];

}


sub transcript_display_xref_sources {
  my $self     = shift;

  my @list = qw(RFAM
	      miRBase
	      RGD 
	      Uniprot/SWISSPROT
	      Uniprot/Varsplic
);


  my %ignore;

  return [\@list,\%ignore];
  
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

package XrefMapper::mus_musculus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };
use strict;

sub get_set_lists {

  return [["ExonerateGappedBest1", ["mus_musculus","*"]]];

}

sub consortium {

  return "MGI";

}

sub gene_description_sources {

  return ("miRBase",
	  "RFAM", 
          "IMGT/GENE_DB",
	  "MGI",
	  "Uniprot/SWISSPROT", 
	  "Uniprot/Varsplic", 
	  "RefSeq_peptide", 
	  "RefSeq_dna", 	  
	  "Uniprot/SPTREMBL" );

}

sub special_filter {

  return ('\(?[0-9A-Z]{10}RIK PROTEIN\)?[ \.]',
	  'RIKEN CDNA [0-9A-Z]{10} GENE',
	  '.*RIKEN FULL-LENGTH ENRICHED LIBRARY.*PRODUCT:',
	  '.*RIKEN FULL-LENGTH ENRICHED LIBRARY.*',
	  '\(*HYPOTHETICAL\s+.*',
	  '^UNKNOWN\s+.*',
	  'CDNA SEQUENCE\s?,? [A-Z]+\d+[ \.;]',
	  'CLONE MGC:\d+[ \.;]',
	  ' MGC:\s*\d+[ \.;]',
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


sub gene_description_filter_regexps {

  return ('\(*HYPOTHETICAL\s+.*',
	  '^UNKNOWN\s+.*',
	  'CDNA SEQUENCE\s?,? [A-Z]+\d+[ \.;]',
	  'CLONE MGC:\d+[ \.;]',
	  ' MGC:\s*\d+[ \.;]',
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

sub get_list_of_sources_for_one_max_per_transcript{
  my $self = shift;

  my @list = qw(MGI);

  return @list;
}

1;

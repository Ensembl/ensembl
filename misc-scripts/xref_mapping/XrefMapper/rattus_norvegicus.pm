package XrefMapper::rattus_norvegicus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["rattus_norvegicus","*"]]];

}


sub transcript_display_xref_sources {

my @list = qw(RFAM
	      miRBase
	      HUGO
	      RGD
	      MGI
	      flybase_symbol
	      Anopheles_symbol
	      Genoscope_annotated_gene
	      Uniprot/SWISSPROT
	      Uniprot/Varsplic
	      RefSeq_peptide
	      RefSeq_dna
	      Uniprot/SPTREMBL
	      EntrezGene);


  my %ignore;
  $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';
  
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

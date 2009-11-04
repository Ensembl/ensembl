package XrefMapper::rattus_norvegicus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["rattus_norvegicus","*"]]];

}



sub transcript_display_xref_sources {
  my $self     = shift;
  my $fullmode = shift;

  my @list = qw(RFAM
	      miRBase
	      RGD
	      MGI
	      flybase_symbol
	      Anopheles_symbol
	      Genoscope_annotated_gene
	      Uniprot/SWISSPROT
	      Uniprot/Varsplic
	      Uniprot/SPTREMBL
	      EntrezGene);


  my %ignore;
  
  if(!$fullmode){
    $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';
  }
  else{
    $ignore{"EntrezGene"} = 'select ox.object_xref_id from object_xref ox, dependent_xref dx, source s1, xref x1, source s2, xref x2 where ox.object_xref_id = dx.object_xref_id and dx.dependent_xref_id = x1.xref_id and x1.source_id = s1.source_id and s1.name = "EntrezGene" and x2.xref_id = dx.master_xref_id and x2.source_id = s2.source_id and (s2.name like "Refseq_dna_predicted" or s2.name like "RefSeq_peptide_predicted") and ox.ox_status = "DUMP_OUT"';
    
  }
  
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

package XrefMapper::danio_rerio;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["danio_rerio","*"]]];

}

sub transcript_display_xref_sources {

  my @list = qw(ZFIN_ID
		MGI
		flybase_symbol
		Anopheles_symbol
		Genoscope_annotated_gene
		Uniprot/SWISSPROT
		RefSeq_peptide
		RefSeq_dna
		Uniprot/SPTREMBL
		EntrezGene);

  my %ignore;
  
  return [\@list,\%ignore];

}

sub gene_description_filter_regexps {

  return ();

}

1;

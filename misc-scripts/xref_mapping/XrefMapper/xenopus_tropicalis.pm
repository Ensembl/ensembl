package XrefMapper::xenopus_tropicalis;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists{

  return [["ExonerateGappedBest_90_perc_id", ["xenopus_tropicalis","*"]]];

}

sub transcript_display_xref_sources {

  return ('Xenopus_Jamboree',
	  'MarkerSymbol',
	  'flybase_symbol',
	  'Anopheles_symbol',
	  'Genoscope_annotated_gene',
	  'Genoscope_predicted_transcript',
	  'Genoscope_predicted_gene',
	  'Uniprot/SWISSPROT',
	  'RefSeq_peptide',
	  'RefSeq_dna',
	  'Uniprot/SPTREMBL',
	  'RefSeq_peptide_predicted',
	  'RefSeq_dna_predicted',
	  'EntrezGene');

}

1;

package XrefMapper::xenopus_tropicalis;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists{

  return [["ExonerateGappedBest_90_perc_id", ["xenopus_tropicalis","*"]]];

}

sub transcript_display_xref_sources {

  my @list = qw(Xenopus_Jamboree
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
  $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';
  
  return [\@list,\%ignore];
}

1;

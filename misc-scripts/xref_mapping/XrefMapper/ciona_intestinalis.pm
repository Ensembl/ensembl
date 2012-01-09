package XrefMapper::ciona_intestinalis;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest5", ["ciona_intestinalis","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["ciona_intestinalis","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["ciona_intestinalis","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["ciona_intestinalis","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["ciona_intestinalis","*"]]];

}

sub gene_description_filter_regexps {

  return ();

}

1;

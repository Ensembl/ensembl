package XrefMapper::monodelphis_domestica;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest5", ["monodelphis_domestica","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["monodelphis_domestica","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["monodelphis_domestica","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["monodelphis_domestica","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["monodelphis_domestica","*"]]];

}

sub gene_description_filter_regexps {

  return ();

}

1;

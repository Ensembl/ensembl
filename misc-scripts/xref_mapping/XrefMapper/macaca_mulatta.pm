package XrefMapper::macaca_mulatta;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest5", ["macaca_mulatta","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["macaca_mulatta","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["macaca_mulatta","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["macaca_mulatta","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["macaca_mulatta","*"]]];

}

sub gene_description_filter_regexps {

  return ();

}

1;

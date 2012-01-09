package XrefMapper::drosophila_erecta;
use strict;

use  XrefMapper::drosophila;
use vars '@ISA';
@ISA = qw{ XrefMapper::drosophila };

sub get_set_lists {

  return [["ExonerateGappedBest5", ["drosophila_erecta","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_erecta","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["drosophila_erecta","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_erecta","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["drosophila_erecta","*"]]];

}

1;

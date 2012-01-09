package XrefMapper::drosophila_virilis;
use strict;

use  XrefMapper::drosophila;
use vars '@ISA';
@ISA = qw{ XrefMapper::drosophila };

sub get_set_lists {

  return [["ExonerateGappedBest5", ["drosophila_virilis","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_virilis","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["drosophila_virilis","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_virilis","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["drosophila_virilis","*"]]];

}

1;

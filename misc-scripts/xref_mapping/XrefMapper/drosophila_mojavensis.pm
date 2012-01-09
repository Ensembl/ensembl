package XrefMapper::drosophila_mojavensis;
use strict;

use  XrefMapper::drosophila;
use vars '@ISA';
@ISA = qw{ XrefMapper::drosophila };


sub get_set_lists {

  return [["ExonerateGappedBest5", ["drosophila_mojavensis","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_mojavensis","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["drosophila_mojavensis","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_mojavensis","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["drosophila_mojavensis","*"]]];

}


1;

package XrefMapper::drosophila_persimilis;
use strict;

use  XrefMapper::drosophila;
use vars '@ISA';
@ISA = qw{ XrefMapper::drosophila };

sub get_set_lists {

  return [["ExonerateGappedBest5", ["drosophila_persimilis","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_persimilis","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["drosophila_persimilis","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_persimilis","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["drosophila_persimilis","*"]]];

}

1;

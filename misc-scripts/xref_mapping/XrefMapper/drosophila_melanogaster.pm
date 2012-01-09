package XrefMapper::drosophila_melanogaster;
use strict;

use  XrefMapper::drosophila;
use vars '@ISA';
@ISA = qw{ XrefMapper::drosophila };

sub get_set_lists {

  return [["ExonerateGappedBest5", ["drosophila_melanogaster","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_melanogaster","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["drosophila_melanogaster","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["drosophila_melanogaster","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["drosophila_melanogaster","*"]]];

}

1;

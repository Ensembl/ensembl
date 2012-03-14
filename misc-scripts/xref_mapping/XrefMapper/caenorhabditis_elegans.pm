package XrefMapper::caenorhabditis_elegans;

use  XrefMapper::wormbase;

use vars qw(@ISA);

@ISA = qw(XrefMapper::wormbase);


sub get_set_lists {

  return [["ExonerateGappedBest5", ["caenorhabditis_elegans","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["caenorhabditis_elegans","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["caenorhabditis_elegans","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["caenorhabditis_elegans","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["caenorhabditis_elegans","*"]]];

}


1;

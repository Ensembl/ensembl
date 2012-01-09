package XrefMapper::canis_familiaris;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest5", ["canis_familiaris","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["canis_familiaris","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["canis_familiaris","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["canis_familiaris","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["canis_familiaris","*"]]];

}

sub gene_description_filter_regexps {

  return ();

}

1;

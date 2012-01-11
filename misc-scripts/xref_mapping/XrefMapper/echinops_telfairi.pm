package XrefMapper::echinops_telfairi;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {
  return [["ExonerateGappedBest_100_perc_id", ["echinops_telfairi","Uniprot/SWISSPROT"]],
	  ["ExonerateGappedBest_100_perc_id", ["echinops_telfairi","Uniprot/SPTREMBL"]],
          ["ExonerateGappedBest1", ["echinops_telfairi","*"]] ];
}

1;

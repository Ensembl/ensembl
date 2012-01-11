package XrefMapper::dasypus_novemcinctus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {
  return [["ExonerateGappedBest_100_perc_id", ["dasypus_novemcinctus","Uniprot/SWISSPROT"]],
	  ["ExonerateGappedBest_100_perc_id", ["dasypus_novemcinctus","Uniprot/SPTREMBL"]],
          ["ExonerateGappedBest1", ["dasypus_novemcinctus","*"]] ];
}

1;

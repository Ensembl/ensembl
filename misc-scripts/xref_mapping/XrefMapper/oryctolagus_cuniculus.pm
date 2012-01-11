package XrefMapper::oryctolagus_cuniculus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {
  return [["ExonerateGappedBest_100_perc_id", ["oryctolagus_cuniculus","Uniprot/SWISSPROT"]],
	  ["ExonerateGappedBest_100_perc_id", ["oryctolagus_cuniculus","Uniprot/SPTREMBL"]],
          ["ExonerateGappedBest1", ["oryctolagus_cuniculus","*"]] ];
}

1;

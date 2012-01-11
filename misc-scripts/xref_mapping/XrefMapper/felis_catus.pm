package XrefMapper::felis_catus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {
  return [["ExonerateGappedBest_100_perc_id", ["felis_catus","Uniprot/SWISSPROT"]],
	  ["ExonerateGappedBest_100_perc_id", ["felis_catus","Uniprot/SPTREMBL"]],
          ["ExonerateGappedBest1", ["felis_catus","*"]] ];
}

1;

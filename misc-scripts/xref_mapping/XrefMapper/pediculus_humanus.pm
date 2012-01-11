package XrefMapper::pediculus_humanus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {
  return [["ExonerateGappedBest_100_perc_id", ["pediculus_humanus","Uniprot/SWISSPROT"]],
	  ["ExonerateGappedBest_100_perc_id", ["pediculus_humanus","Uniprot/SPTREMBL"]],
          ["ExonerateGappedBest1", ["pediculus_humanus","*"]] ];
}

1;

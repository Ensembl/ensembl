package XrefMapper::xenopus_tropicalis;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists{

  return [["ExonerateGappedBest_90_perc_id", ["xenopus_tropicalis","*"]]];

}

1;

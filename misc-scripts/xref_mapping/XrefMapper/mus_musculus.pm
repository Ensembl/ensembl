package XrefMapper::mus_musculus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists{

  return [["ExonerateGappedBest1", ["mus_musculus","*"]]];

}

1;

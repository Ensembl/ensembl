package XrefMapper::homo_sapiens;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists{

  return [["ExonerateGappedBest1", ["homo_sapiens","*"]]];

}

1;

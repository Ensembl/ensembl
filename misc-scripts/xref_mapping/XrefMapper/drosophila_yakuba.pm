package XrefMapper::drosophila_yakuba;
use strict;
use  XrefMapper::drosophila;

use vars '@ISA';
@ISA = qw{ XrefMapper:: XrefMapper::drosophila };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["drosophila_yakuba","*"]]];

}

1;

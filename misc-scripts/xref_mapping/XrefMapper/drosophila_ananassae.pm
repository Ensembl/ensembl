package XrefMapper::drosophila_ananassae;
use strict;

use  XrefMapper::drosophila;
use vars '@ISA';
@ISA = qw{ XrefMapper::drosophila };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["drosophila_ananassae","*"]]];

}


1;

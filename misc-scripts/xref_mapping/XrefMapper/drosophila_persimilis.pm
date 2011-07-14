package XrefMapper::drosophila_persimilis;
use strict;

use  XrefMapper::drosophila;
use vars '@ISA';
@ISA = qw{ XrefMapper::drosophila };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["drosophila_persimilis","*"]]];

}

1;

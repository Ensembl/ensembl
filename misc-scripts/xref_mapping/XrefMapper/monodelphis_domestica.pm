package XrefMapper::monodelphis_domestica;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["monodelphis_domestica","*"]]];

}

sub gene_description_filter_regexps {

  return ();

}

1;

package XrefMapper::drosophila_melanogaster;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["drosophila_melanogaster","*"]]];

}

sub gene_description_filter_regexps {

  return ();

}

1;

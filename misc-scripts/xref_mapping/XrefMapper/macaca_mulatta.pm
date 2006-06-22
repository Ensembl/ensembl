package XrefMapper::macaca_mulatta;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["macaca_mulatta","*"]]];

}

sub gene_description_filter_regexps {

  return ();

}

1;

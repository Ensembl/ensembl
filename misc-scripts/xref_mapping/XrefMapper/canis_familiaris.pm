package XrefMapper::canis_familiaris;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["canis_familiaris","*"]]];

}

sub gene_description_filter_regexps {

  return ();

}

1;

package XrefMapper::fugu_rubripes;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["fugu_rubripes","*"]]];

}

sub gene_description_filter_regexps {

  return ('^HYPOTHETICAL\s+PROTEIN\.?',
	  '^\s*\(?FRAGMENT\)?\.?\s*',
	  '^SIMILAR TO HUMAN CDNA KIAA\d+\.?');

}

1;

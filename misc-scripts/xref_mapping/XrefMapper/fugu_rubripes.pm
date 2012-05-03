package XrefMapper::fugu_rubripes;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };


sub gene_description_filter_regexps {

  return ('^HYPOTHETICAL\s+PROTEIN\.?',
	  '^\s*\(?FRAGMENT\)?\.?\s*',
	  '^SIMILAR TO HUMAN CDNA KIAA\d+\.?');

}

1;

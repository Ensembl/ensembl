package XrefMapper::caenorhabditis_elegans;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["caenorhabditis_elegans","*"]]];

}

sub gene_description_filter_regexps {

  return ('[0-9A-Z]+\.\d*[A-Z]* PROTEIN[ \.]',
	  '\(\d[A-Z]\d+\)\.',
	  '\([0-9A-Z]+\.\d*[A-Z]* PROTEIN\)[ \.]',
	  '^\(*HYPOTHETICAL\s+.*',
	  '^\s*\(FRAGMENT\)\.?\s*$' );

}

1;

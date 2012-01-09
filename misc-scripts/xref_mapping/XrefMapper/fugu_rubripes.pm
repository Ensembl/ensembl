package XrefMapper::fugu_rubripes;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest5", ["fugu_rubripes","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["fugu_rubripes","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["fugu_rubripes","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["fugu_rubripes","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["fugu_rubripes","*"]]];

}

sub gene_description_filter_regexps {

  return ('^HYPOTHETICAL\s+PROTEIN\.?',
	  '^\s*\(?FRAGMENT\)?\.?\s*',
	  '^SIMILAR TO HUMAN CDNA KIAA\d+\.?');

}

1;

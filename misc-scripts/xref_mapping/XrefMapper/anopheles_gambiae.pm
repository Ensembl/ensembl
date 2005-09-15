package XrefMapper::homo_sapiens;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["anopheles_gambiae","*"]]];

}

# transcript, gene display_xrefs and descriptions can use defaults
# since anopheles_symbol is "before" Uniprot

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {

  return ();

}


1;

package XrefMapper::anopheles_gambiae;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["anopheles_gambiae","*"]]];

}

# transcript, gene display_xrefs can use defaults
# since anopheles_symbol is "before" Uniprot

# If there is an Anopheles_symbol xref, use its description
sub gene_description_sources {

  return ("RefSeq_dna_predicted",
	  "RefSeq_peptide_predicted",
	  "Uniprot/SPTREMBL",
	  "RefSeq_dna",
	  "RefSeq_peptide",
	  "Uniprot/SWISSPROT",
	  "Anopheles_symbol");

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {

  return ();

}


1;

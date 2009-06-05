

package XrefMapper::anopheles_gambiae;

use  XrefMapper::BasicMapper;
use  XrefMapper::VBCoordinateMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1_agam", ["anopheles_gambiae","*"]]];

}

# transcript, gene display_xrefs can use defaults
# since anopheles_symbol is "before" Uniprot

# If there is an Anopheles_symbol xref, use its description

# mh4 says Anopheles_symbol doesn't get chosen over UniP
# (but they do get chosen in other cases)

sub gene_description_sources {

  return ("VB_Community_Symbol",
	  "Uniprot/SWISSPROT",
	  #"RefSeq_peptide",
	  #"RefSeq_dna",
	  #"Uniprot/SPTREMBL",
	  #"RefSeq_peptide_predicted",
	  #"RefSeq_dna_predicted",
	  #"EntrezGene");
          );
}

sub transcript_display_xref_sources {

  my @list = qw(RFAM
		miRBase
		VB_Community_Symbol
		Uniprot/SWISSPROT
		#Uniprot/Varsplic
		#Uniprot/SPTREMBL);
   );

  my %ignore;
  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {

  return ();

}

1;

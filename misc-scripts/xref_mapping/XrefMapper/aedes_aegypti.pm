package XrefMapper::aedes_aegypti;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {
  return [["ExonerateGappedBest1", ["aedes_aegypti","*"]]];
}

#Reverse order: last one has higher priority!
sub gene_description_sources {
  return ("VB_External_Description",
	  "Uniprot/SWISSPROT",
          "VB_Community_Annotation"
	  );
}

sub transcript_display_xref_sources {
  my @list = qw(
		Uniprot/SWISSPROT
		VB_Community_Annotation
                );

  my %ignore;
  #$ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';

  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {
  return ();
}1;

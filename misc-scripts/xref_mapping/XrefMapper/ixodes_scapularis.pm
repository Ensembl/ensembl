package XrefMapper::ixodes_scapularis;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {
  return [["ExonerateGappedBest1_agam", ["ixodes_scapularis","*"]]];  #use same parameters as Anoph. and Culex
}

#Reverse order: second one has higher priority!
sub gene_description_sources {
  return ("Uniprot/SWISSPROT",
	  "Ixodes_ManualAnnotation",
  ) ;
}

#Reverse order: second one has higher priority!
sub transcript_display_xref_sources {
  my @list = qw(Uniprot/SWISSPROT
		Ixodes_ManualAnnotation
                );

  my %ignore;
  $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';

  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {
  return ();
}1;

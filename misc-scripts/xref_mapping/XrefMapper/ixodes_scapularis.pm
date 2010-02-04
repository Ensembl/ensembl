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
  my $self     = shift;
  my $fullmode = shift;

  my @list = qw(Uniprot/SWISSPROT
		Ixodes_ManualAnnotation
                );

  my %ignore;
  if(!$fullmode){
    $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';
  }
  else{
    $ignore{"EntrezGene"} = 'select ox.object_xref_id from object_xref ox, dependent_xref dx, source s1, xref x1, source s2, xref x2 where ox.object_xref_id = dx.object_xref_id and dx.dependent_xref_id = x1.xref_id and x1.source_id = s1.source_id and s1.name = "EntrezGene" and x2.xref_id = dx.master_xref_id and x2.source_id = s2.source_id and (s2.name like "Refseq_dna_predicted" or s2.name like "RefSeq_peptide_predicted") and ox.ox_status = "DUMP_OUT"';

  }

  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {
  return ();
}1;

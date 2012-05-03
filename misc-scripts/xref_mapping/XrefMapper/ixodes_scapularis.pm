package XrefMapper::ixodes_scapularis;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };


sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1_55_perc_id';
  my %override_method_for_source = ( );

  return $default_method, \%override_method_for_source;
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

  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {
  return ();
}1;

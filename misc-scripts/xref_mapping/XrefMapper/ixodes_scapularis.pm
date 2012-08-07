package XrefMapper::ixodes_scapularis;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };


sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1_55_perc_id';
  my %override_method_for_source = ( );

  return $default_method, \%override_method_for_source;
}


sub gene_description_sources {
  return ("Ixodes_ManualAnnotation",
	  "Uniprot/SWISSPROT",	  
  ) ;
}

sub gene_display_xref_sources {
  my $self     = shift;
  my $fullmode = shift;

  my @list = qw(Ixodes_ManualAnnotation
                Uniprot/SWISSPROT		
                );

  my %ignore;

  return [\@list,\%ignore];

}

sub transcript_display_xref_sources {
  my $self     = shift;
  my $fullmode = shift;

  my @list = qw(Ixodes_ManualAnnotation
                Uniprot/SWISSPROT		
                );

  my %ignore;

  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {
  return ();
}1;

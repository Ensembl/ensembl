package XrefParser::Vega_TranParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of Vega transcript name records and assign direct xrefs
# All assumed to be linked to transcripts

sub run {

  my ($self, $source_id, $species_id, $file) = @_;

  die "No longer used try HGNC_curated_transcript\n";

  my $vega_io = $self->get_filehandle($file);

  my $clone_source_id =
    $self->get_source_id_for_source_name('Clone_based_vega_transcript');
  my $curated_source_id =
    $self->get_source_id_for_source_name('HGNC_curated_transcript');
  
  if ( !defined $vega_io ) {
    print STDERR "Could not open $file\n";
    return 1;
  }

  my $xref_count = 0;
  my $line_count = 0;
  
  my %vega_name;
  my %vega_priority;

  while ( $_ = $vega_io->getline() ) {
    my $line_count++;
    my ($stable_id, $name, $priority) = split;
    
    if(!defined($vega_name{$stable_id}) or $vega_priority{$stable_id} < $priority){
      $vega_name{$stable_id} = $name;
      $vega_priority{$stable_id} = $priority;
    }	
  }
  

  foreach my $stable_id (keys %vega_name){
    my $name = $vega_name{$stable_id};
    my $id = $curated_source_id;
    if($name =~ /[.]/){
      $id = $clone_source_id;
    }
    my $xref_id = $self->add_xref($name, "" , $name , "", $id, $species_id);
    $xref_count++;
    
    
    $self->add_direct_xref($xref_id, $stable_id, "transcript", "");
  }

  print "Parsed $line_count lines from $file, added $xref_count xrefs and $xref_count direct_xrefs\n" if($verbose);

  $vega_io->close();

  return 0;
}

1;

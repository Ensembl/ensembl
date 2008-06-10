package XrefParser::Vega_TranParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of Vega transcript name records and assign direct xrefs
# All assumed to be linked to transcripts

sub run {

  my ($self, $source_id, $species_id, $file) = @_;

  my $vega_io = $self->get_filehandle($file);

  if ( !defined $vega_io ) {
    print "Could not open $file\n";
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
    my $xref_id = $self->add_xref($name, "" , $name , "", $source_id, $species_id);
    $xref_count++;
    
    
    $self->add_direct_xref($xref_id, $stable_id, "transcript", "");
  }

  print "Parsed $line_count lines from $file, added $xref_count xrefs and $xref_count direct_xrefs\n";

  $vega_io->close();

  return 0;
}

1;

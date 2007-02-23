package XrefParser::OTTTParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of Ensembl - Vega OTTT transcript mappings
# ENST00000373795:	OTTHUMT00000010392
# ENST00000374603:	OTTHUMT00000057024
# ENST00000372604:	OTTHUMT00000057746
# ENST00000329151:	OTTHUMT00000011475

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my $ottt_io = $self->get_filehandle($file);

  if ( !defined $ottt_io ) {
    print "Could not open $file\n";
    return 1;
  }

  my $line_count = 0;
  my $xref_count = 0;

  my $xref_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$source_id AND species_id=$species_id");

  while ( $_ = $ottt_io->getline() ) {
    my ($ens, $ottt) = split;

    $ens =~ s/://g;

    $line_count++;

    # check if an xref already exists
    $xref_sth->execute($ottt);
    my $xref_id = ($xref_sth->fetchrow_array())[0];
    if (!$xref_id) {
      $xref_id = $self->add_xref($ottt, '', $ottt, "", $source_id, $species_id);
      $xref_count++;
    }

    $self->add_direct_xref($xref_id, $ens, "transcript", "");

  }

  $ottt_io->close();

  print "Parsed $line_count OTTT identifiers from $file, added $xref_count xrefs and $line_count direct_xrefs\n";


  return 0;
}

1;

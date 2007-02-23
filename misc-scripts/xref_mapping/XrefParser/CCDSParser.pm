package XrefParser::CCDSParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of CCDS records and assign direct xrefs
# All assumed to be linked to transcripts
# The same CCDS may be linked to more than one transcript, but need to only
# add the xref once, so check if it already exists before adding it.

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my $ccds_io = $self->get_filehandle($file);

  if ( !defined $ccds_io ) {
      print "Could not open $file\n";
      return 1;
  }

  my $line_count = 0;
  my $xref_count = 0;

  my $xref_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND version=? AND source_id=$source_id AND species_id=$species_id");

  while ( $_ = $ccds_io->getline() ) {
    my ($stable_id, $ccds) = split;

    my ($acc, $version) = split (/\./, $ccds);
    $line_count++;

    # check if an xref already exists
    $xref_sth->execute($acc, $version);
    my $xref_id = ($xref_sth->fetchrow_array())[0];
    if (!$xref_id) {
      $xref_id = $self->add_xref($acc, $version, $ccds, "", $source_id, $species_id);
      $xref_count++;
    }

    $self->add_direct_xref($xref_id, $stable_id, "transcript", "");

  }

  print "Parsed $line_count CCDS identifiers from $file, added $xref_count xrefs and $line_count direct_xrefs\n";

  $ccds_io->close();
  return 0;
}

1;

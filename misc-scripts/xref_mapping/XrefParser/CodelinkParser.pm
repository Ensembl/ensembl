package XrefParser::CodelinkParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

# Parser for Codelink probes

#>GE469530
#TTGTTTTCAGCTTGCTTCTGTCATTCTTCC
#>GE469548
#CACAGTTGGGTGAAGCTGGTGATGAAGGTA

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my @xrefs;

  local $/ = "\n>";

  my $codelink_io = $self->get_filehandle($file);
  if ( !defined $codelink_io ) {
    print "ERROR: Could not open $file\n";
    return 1;    # 1 = error
  }

  while ( $_ = $codelink_io->getline() ) {
    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header - only accession for now
    my $accession = $header;

    # make sequence into one long string - probably not necessary for short probes
    $sequence =~ s/\n//g;

    # build the xref object and store it
    $xref->{ACCESSION}     = $accession;
    $xref->{LABEL}         = $accession;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = 'dna';
    $xref->{STATUS}        = 'experimental';

    push @xrefs, $xref;

  }

  $codelink_io->close();

  print scalar(@xrefs) . " Codelink xrefs succesfully parsed\n";

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print "Done\n";
  return 0; #successful
}

1;

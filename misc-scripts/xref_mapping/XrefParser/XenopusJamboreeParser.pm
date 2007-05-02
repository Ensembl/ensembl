package XrefParser::XenopusJamboreeParser;

# Parse annotated peptides from Xenopus Jamboree

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

# Xenopus Jamboree peptides file format: fasta, e.g.

# >anxa2 BC061610
# MSLIHEILGKLSLEGNQSSTRQSTLGSVKASSNFDAERDAAAIETAIKTKGVDELTII


sub run {

  my ($self, $source_id, $species_id, $file) = @_;

  my @xrefs;

  local $/ = "\n>";

  my $file_io = $self->getline($file);

  if ( !defined $file_io ) {
    print "ERROR: Could not open $file\n";
    return 1;    # 1 error
  }

  while ( $_ = $file_io->getline() ) {
    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header
    my ($accession, @other) = split /\s+/, $header;
    my $other = join(" ", @other);

    # make sequence into one long string
    $sequence =~ s/\n//g;

    # build the xref object and store it
    $xref->{ACCESSION}     = $accession;
    $xref->{LABEL}         = $accession;
    $xref->{DESCRIPTION}   = '';
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS}        = 'experimental';

    push @xrefs, $xref;

  }

  $file_io->close();

  print scalar(@xrefs) . " XenopusJamboreeParser xrefs succesfully parsed\n";

  if(!defined(XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs))){
    return 1; #1 error
  }

  print "Done\n";
  return 0;
}

1;

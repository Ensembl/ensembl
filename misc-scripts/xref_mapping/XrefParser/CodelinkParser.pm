package XrefParser::CodelinkParser;

use strict;
use File::Basename;
use Carp;
use base qw( XrefParser::BaseParser );

# Parser for Codelink probes

#>GE469530
#TTGTTTTCAGCTTGCTTCTGTCATTCTTCC
#>GE469548
#CACAGTTGGGTGAAGCTGGTGATGAAGGTA

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) or (!defined $release_file)){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my @xrefs;

  local $/ = "\n>";

  my $codelink_io = $self->get_filehandle($file);
  if ( !defined $codelink_io ) {
    print STDERR "ERROR: Could not open $file\n";
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


  $self->upload_xref_object_graphs(\@xrefs);

  print scalar(@xrefs) . " Codelink xrefs succesfully parsed\n" if($verbose);

  return 0; #successful
}

1;

package XrefParser::IlluminaParser;

use strict;

use base qw( XrefParser::BaseParser );

# Parser for Illumina V2 xrefs - V1 are done by the vanilla FastaParser

# Search_key,Target,ProbeId,Gid,Transcript,Accession,Symbol,Type,Start,Probe_Sequence,Definition,Ontology,Synonym
# ILMN_89282,ILMN_89282,0004760445,23525203,Hs.388528,BU678343,"",S,349,CTCTCTAAAGGGACAACAGAGTGGACAGTCAAGGAACTCCACATATTCAT,"UI-CF-EC0-abi-c-12-0-UI.s1 UI-CF-EC0 Homo sapiens cDNA clone UI-CF-EC0-abi-c-12-0-UI 3, mRNA sequence",,
# ILMN_35826,ILMN_35826,0002760372,89042416,XM_497527.2,XM_497527.2,"LOC441782",S,902,GGGGTCAAGCCCAGGTGAAATGTGGATTGGAAAAGTGCTTCCCTTGCCCC,"PREDICTED: Homo sapiens similar to spectrin domain with coiled-coils 1 (LOC441782), mRNA.",,

# Note that "definition" column often has commas.

sub run {

  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  my @xrefs;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "Could not open $file\n";
    return 1;
  }

  while ( $_ = $file_io->getline() ) {
    chomp;

    my $xref;

    # strip ^M at end of line
    $_ =~ s/\015//g;

    my @bits = split(/,[^ ]/);
    my $illumina_id = $bits[0];
    next if ($illumina_id eq "Search_key");   # skip header
    next if (!$illumina_id); # skip lines with missing accessions

    my $sequence = $bits[9];

    my $type = $bits[7];
    # XXX what about "type" column?


    my ($description) = $bits[10];
    $description =~ s/\"//g;

    # build the xref object and store it
    $xref->{ACCESSION}     = $illumina_id;
    $xref->{LABEL}         = $illumina_id;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{DESCRIPTION}   = $description;
    $xref->{SEQUENCE_TYPE} = 'dna';
    $xref->{STATUS}        = 'experimental';

    push @xrefs, $xref;

  }

  $file_io->close();

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print scalar(@xrefs) . " Illumina V2 xrefs succesfully parsed\n" if($verbose);


  return 0;
}

1;

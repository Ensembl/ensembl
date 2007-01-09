package XrefParser::IlluminaParser;

use strict;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# Parser for Illumina V2 xrefs - V1 are done by the vanilla FastaParser

# Search_key,Target,ProbeId,Gid,Transcript,Accession,Symbol,Type,Start,Probe_Sequence,Definition,Ontology,Synonym
# ILMN_89282,ILMN_89282,0004760445,23525203,Hs.388528,BU678343,"",S,349,CTCTCTAAAGGGACAACAGAGTGGACAGTCAAGGAACTCCACATATTCAT,"UI-CF-EC0-abi-c-12-0-UI.s1 UI-CF-EC0 Homo sapiens cDNA clone UI-CF-EC0-abi-c-12-0-UI 3, mRNA sequence",,
# ILMN_35826,ILMN_35826,0002760372,89042416,XM_497527.2,XM_497527.2,"LOC441782",S,902,GGGGTCAAGCCCAGGTGAAATGTGGATTGGAAAAGTGCTTCCCTTGCCCC,"PREDICTED: Homo sapiens similar to spectrin domain with coiled-coils 1 (LOC441782), mRNA.",,

# Note that "definition" column often has commas.

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my @xrefs;

  if(!open(FILE,"<".$file)){
    print "Could not open $file\n";
    return 1;
  }

  while (<FILE>) {

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

  close(FILE);

  print scalar(@xrefs) . " Illumina V2 xrefs succesfully parsed\n";

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);


  return 0;
}


sub new {

  my $self = {};
  bless $self, "XrefParser::IlluminaParser";
  return $self;

}

1;

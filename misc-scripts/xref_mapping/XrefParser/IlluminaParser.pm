package XrefParser::IlluminaParser;

use strict;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# Parser for Illumina refs

#Gid,Accession,Symbol,Probe_Sequence,Definition
#GI_4505876,NM_000445.1,PLEC1,AACACTAACCTGACCGTGGGCGGGGCCTTGCGGTATCCGCCCCCAATAAA,"Homo sapiens plectin 1, intermediate filament binding protein 500kDa (PLEC1), transcript variant 1, mRNA."
#GI_18860893,NM_130393.1,PTPRD,CTACAGGCCCTTCAATATCCATGGAGTCTCTTCTGAGCCATACAGGGCAC,"Homo sapiens protein tyrosine phosphatase, receptor type, D (PTPRD), transcript variant 4, mRNA."

# Note that "defintion" column often has commas.

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

    my ($accession, $refseq, $symbol, $sequence) = split(/,/); # ignore description for now
    next if ($accession eq "Gid" && $refseq eq "Accession");   # skip header

    next if (!$accession); # skip lines with missing accessions and "null" for RefSeq

    my ($description) = $_ =~ /.*\"([^\"]*)\.\"/; # parse description

    # build the xref object and store it
    $xref->{ACCESSION}     = $accession;
    $xref->{LABEL}         = $accession;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{DESCRIPTION}   = $description;
    $xref->{SEQUENCE_TYPE} = 'dna';
    $xref->{STATUS}        = 'experimental';

    push @xrefs, $xref;

  }

  close(FILE);

  print scalar(@xrefs) . " Illumina xrefs succesfully parsed\n";

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);


  return 0;
}


sub new {

  my $self = {};
  bless $self, "XrefParser::IlluminaParser";
  return $self;

}

1;

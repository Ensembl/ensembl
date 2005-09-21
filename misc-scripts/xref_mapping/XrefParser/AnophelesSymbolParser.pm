package XrefParser::AnophelesSymbolParser;

use strict;
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# AnophelesSymbol database dump for anopheles - FASTA format
#
# >ANXB10B,ANXB10B,ANXB10B,protein,Annexin B10B
# MSWYYTPHPTVVPAEDFDASADANALRKAMKGFGTDEQAIIDILCARSNGQRQEIAEAFKRELGRDLIDDLKSELGGKFEDVILGLMLRPEAYLCKQLHKAMDGIGTDEKSLIEII
# CPQTNDQIRAIVDCYEEMYSRPLAEHLCSETSGSFRRLLTMIIVGSRDPQGTVDPELAVEQAKQLYDAGEGKLGTDEEVFYKILAHASFDQLEIVFEEYKSLSGRTIEQALKAELS
# GELYDALSAIVECVQMAPHFFAKRLHKAMDGVGTDDATLIRIIVSRSEIDLQNIKDEFEQMYNKTLVSAVRSETSGDYKRALCALIGNA

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  next if (/^File:/);   # skip header

  my @xrefs;

  local $/ = "\n>";

  open(FILE,"<".$file) || die "Could not open $file\n";

  while (<FILE>) {

    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header - just use first part
    my ($accession, @rest) = split /,/, $header;

    # make sequence into one long string
    $sequence =~ s/\n//g;

    # build the xref object and store it
    $xref->{ACCESSION}     = $accession;
    $xref->{LABEL}         = $accession;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS}        = 'experimental';

    push @xrefs, $xref;

  }

  close (FILE);

  print scalar(@xrefs) . " AnophelesSymbol xrefs succesfully parsed\n";

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print "Done\n";

}


sub new {

  my $self = {};
  bless $self, "XrefParser::AnophelesSymbolParser";
  return $self;

}

1;

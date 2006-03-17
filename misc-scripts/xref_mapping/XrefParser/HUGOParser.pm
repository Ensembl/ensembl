package XrefParser::HUGOParser;

use strict;
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);
my $xref_sth ;
my $dep_sth;
my $syn_sth;

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: HUGOParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run(@ARGV);
}

sub run {

  my $self = shift if (defined(caller(1)));
  my $file = shift;

  my $source_id = shift;
  my $species_id = shift;

  print STDERR "source = $source_id\tspecies = $species_id\n";
  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }


  my (%swiss)  =  %{XrefParser::BaseParser->get_valid_codes("uniprot",$species_id)};
  my (%refseq) =  %{XrefParser::BaseParser->get_valid_codes("refseq",$species_id)};

  my $swiss_count = 0;
  my $refseq_count = 0;
  my $mismatch = 0;

  if(!open (HUGO, "<$file")){
    print  "ERROR: Can't open HUGO file $file\n";
    return 1;
  }

  <HUGO>;

  #23	ABAT	4-aminobutyrate aminotransferase		P80404
  #29	ABCA1	ATP-binding cassette, sub-family A (ABC1), member 1	ABC1, HDLDT1	O95477
  #40	ABCB1	ATP-binding cassette, sub-family B (MDR/TAP), member 1	PGY1, MDR1, CLCS	P-gp, CD243, GP170, ABC20	P08183	NM_000927

  while (<HUGO>) {

    chomp;

    # 0 HGNC ID	           # primary accession
    # 1 Approved Symbol    # label
    # 2 Approved Name      # description
    # 3 Previous Symbols   # synonyms
    # 4 Aliases            # aliases
    # 5 UniProt ID         # uniprot accession
    # 6 RefSeq ID

    my @array = split(/\t/,$_);

    # Use the RefSeq if available as this is manually curated
    # If no RefSeq, use the Swissprot instead

    my $master;
    if ($array[6]) {               # RefSeq
      $master = $refseq{$array[6]};
      $refseq_count++;
    } elsif ($array[5]) {        # Uniprot
      $master = $swiss{$array[5]};
      $swiss_count++;
    }

    if (!$master) {
      $mismatch++;
      next;
    }

    #print $array[5] ." " . $array[6] . " " . $master . " " . $swiss_count . " " . $refseq_count . " " . $mismatch . "\n";

    XrefParser::BaseParser->add_to_xrefs($master, $array[0], '', $array[1], $array[2], "", $source_id, $species_id);

    if (defined($array[3])) {     # dead name, add to synonym
      my @array2 = split(',\s*', $array[3]);
      foreach my $arr (@array2){
	XrefParser::BaseParser->add_to_syn($array[0], $source_id, $arr);
      }
    }

    if (defined($array[4])) {     # alias, add to synonym
      my @array2 = split(',\s*', $array[4]);
      foreach my $arr (@array2){
	XrefParser::BaseParser->add_to_syn($array[0], $source_id, $arr);
      }
    }

  } # while HUGO

  close (HUGO);

  print "Loaded a total of " . ($swiss_count + $refseq_count) . " HUGO xrefs, $refseq_count from RefSeq curated mappings and $swiss_count from Uniprot (mapped)\n";

  print "$mismatch xrefs could not be associated via RefSeq or Uniprot\n";

  return 0; # successful

}

sub rename_url_file{
  return "hugo.txt";
}

sub new {

  my $self = {};
  bless $self, "XrefParser::HUGOParser";
  return $self;

}

1;
    


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
#  my (%refseq) =  %{XrefParser::BaseParser->get_valid_codes("refseq",$species_id)};

  my $count = 0;
  my $mismatch = 0;
  open (HUGO, "<$file") || die "Can't open hugo file $file\n";


  <HUGO>;
#23	ABAT	4-aminobutyrate aminotransferase		P80404
#29	ABCA1	ATP-binding cassette, sub-family A (ABC1), member 1	ABC1, HDLDT1	O95477
  while (<HUGO>) {
    chomp;
# 0 HGNC ID	       # primary accession
# 1 Approved Symbol    # label
# 2 Approved Name      # description
# 3 Previous Symbols   # synonyms
# 4 UniProt ID         # uniprot accession
    my @array = split(/\t/,$_);

    if ($array[4]) {  #swissprot
      my $master = $swiss{$array[4]};
      if(!defined($master)){
	$mismatch++;
      }
      else{
	XrefParser::BaseParser->add_to_xrefs($master,$array[0],'',$array[1],$array[2],"",$source_id,$species_id);
	$count++;
	if(defined($array[3])){ #dead name add to synonym
	  my @array2 = split(',\s*',$array[3]);
	  foreach my $arr (@array2){
#	    print "adding synonym ".$arr." for ".$hugo{$hgnc}." ($hgnc)\n";
	    XrefParser::BaseParser->add_to_syn($array[0], $source_id, $arr);
	  }
	}
      }
      #	print "$array[1]\tSPTR\t$hgnc\tHUGO\t$hugo_id{$hgnc}\t$hugo_syn{$hgnc}\tXREF\n";
    }

  }
  close (HUGO);
  print "\t$count xrefs succesfully loaded\n";
  print "\t$mismatch xrefs ignored\n";
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
    


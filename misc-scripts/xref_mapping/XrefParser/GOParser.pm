package XrefParser::GOParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);



# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: GoParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0]);

}

sub run {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id = shift;
  my $species_id = shift;

  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }


  my (%swiss) = XrefParser::BaseParser->get_valid_codes("uniprot",$species_id);
  my (%refseq) = XrefParser::BaseParser->get_valid_codes("refseq",$species_id);

  my $count  = 0;

  open(GO,"<".$file) || die "Could not open $file\n";

  my $taxon_line = "taxon:".$species_id;

  while (<GO>) {
    if(/$taxon_line/){
      chomp;
      my @array = split (/\t/,$_);
      $array[9] =~ s/\'/\\\'/g;
      my $master=0;
      if($array[0] =~ /ENSEMBL/){
	#these might be good for a check
	# match GO to Uniprot
	# match Uniprot to ENSEMBL
	# check ENSEMBL's are the same.
      }
      elsif($array[0] =~ /RefSeq/){
	if($refseq{$array[1]}){
	  XrefParser::BaseParser->add_to_xrefs($refseq{$array[1]},$array[4],'',$array[4],'',$array[6],$source_id,$species_id);
	  $count++;
	}
      }
      elsif($array[0] =~ /UniProt/){
	if($swiss{$array[1]}){
	  XrefParser::BaseParser->add_to_xrefs($swiss{$array[1]},$array[4],'',$array[4],'',$array[6],$source_id,$species_id);
	  $count++;
	}
      }
      else{
	print STDERR "unknown type ".$array[0]."\n";
      }
    }
  }
  print "\t$count GO dependent xrefs added\n"; 
}

sub new {

  my $self = {};
  bless $self, "XrefParser::GOParser";
  return $self;

}
 
1;

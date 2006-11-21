package XrefParser::RGDParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

my $xref_sth ;
my $dep_sth;


# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: RGDParser.pm file <source_id> <species_id>\n\n";
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


#  my (%genbank) = %{XrefParser::BaseParser->get_valid_codes("EMBL",$species_id)};
  my (%refseq) = %{XrefParser::BaseParser->get_valid_codes("refseq",$species_id)};
  

  print STDERR "refseq->". %refseq."\n";

  if(!open(RGD,"<".$file)){
    print "ERROR: Could not open $file\n";
    return 1;
  }
  my $line = <RGD>;
  chomp $line;
  my @linearr = split(/\t/,$line);

  #
  #warn if sanity check fails
  #

  if($linearr[0] =~ /GENE_RDB_ID/){
   die ($linearr[0]."!= GENE_RDB_ID is not the first element in the header\n");
  }
  if($linearr[1] ne "SYMBOL"){
    die ("SYMBOL is not the second element in the header\n");
  }
  if($linearr[2] ne "NAME"){
    die ("NAME is not the third element in the header\n");
  }
  if($linearr[16] ne "GENBANK_NUCLEOTIDE"){
    die ("GENBANK_NUCLEOTIDE is not the seventeenth element in the header but ".$linearr[16]." is.\\n");
  }
  
  my $count= 0;
  my $mismatch = 0;
  while ($line = <RGD>) {
    chomp $line;
    my ($rgd, $symbol, $name, $refseq) = (split (/\t/,$line))[0,1,2,16];
    my @nucs = split(/\,/,$refseq);
    my $done = 0;
    my $failed_list ="";
    foreach my $nuc (reverse @nucs){
      if(!$done){
	my $xref=undef; 
	$xref=$refseq{$nuc} if defined($refseq{$nuc});
	if(defined($xref)){
	  $done = 1;
	  XrefParser::BaseParser->add_to_xrefs($xref,"RGD:".$rgd,"",$symbol,$name,"",$source_id,$species_id);
	  $count++;
	}
	else{
	  $failed_list .= " $nuc";
	}
      }
    }
    if(!$done){
      print STDERR "$rgd FAILED for $failed_list\n";
      $self->add_xref("RGD:".$rgd,"",$symbol,$name,$source_id,$species_id);
      $mismatch++;
    }

  }
  print "\t$count xrefs succesfully loaded and dependent on refseq\n";
  print "\t$mismatch xrefs added but with NO dependencies\n";
  return 0;
}


sub new {

  my $self = {};
  bless $self, "XrefParser::RGDParser";
  return $self;

}
 
1;

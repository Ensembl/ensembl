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
    print "\nUsage: RGDParser.pm file\n\n";
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
    print "source id is $source_id \n";
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
    print "species id is $species_id \n";
  }


#  my (%swiss) = XrefParser::BaseParser->get_valid_codes("uniprot",$species_id);
  my (%genbank) = XrefParser::BaseParser->get_valid_codes("EMBL",$species_id);
  my (%refseq) = XrefParser::BaseParser->get_valid_codes("refseq",$species_id);

  open(RGD,"<".$file) || die "Could not open $file\n";
#Genbank_Nucleotide      Gene_Symbol     RGD_ID
#U40064  Ppard   3370
#AF218575        Nbn     621420

  my $count= 0;
  my $mismatch = 0;
  <RGD>;
  while (<RGD>) {
    chomp;
    my @array = split (/\t/,$_);
    my $xref=undef; 
    $xref=$refseq{$array[0]} if defined($refseq{$array[0]});
    $xref=$genbank{$array[0]}  if defined($genbank{$array[0]});
    if(defined($xref)){
      XrefParser::BaseParser->add_to_xrefs($xref,"RGD:".$array[2],"",$array[1],"",$source_id,$species_id);
      $count++;
    }
    else{
      $mismatch++;
    }
  }
  print "$count xrefs succesfully loaded\n";
  print "$mismatch xrefs ignored\n";
}


sub new {

  my $self = {};
  bless $self, "XrefParser::RGDParser";
  return $self;

}
 
1;

package XrefParser::HUGOParser;
  
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

  my $dir = dirname($file);
 
  my %hugo;
    
  open (ENS4, $dir."/ens4.txt") || die "Can't open hugo ens4 $dir/ens4.txt\n";
  #HGNC    Symbol  Literature Aliases      Withdrawn Symbols
  #5       A1BG
  #7       A2M
  <ENS4>;  #header line
  while(<ENS4>){
    chomp;
    my @array = split(/\t/,$_);
    my $hgnc = $array[0];
    my $label = $array[1];
    
    $hugo{$hgnc} = $label
  }
  close ENS4;

  
  my (%swiss)  =  XrefParser::BaseParser->get_valid_codes("uniprot",$species_id);
  my (%refseq) =  XrefParser::BaseParser->get_valid_codes("refseq",$species_id);
 

  my $count = 0;
  my $mismatch = 0;
  open (ENS1, $dir."/ens1.txt") || die "Can't open hugo ens1  $dir/ens1.txt\n";
  #HGNC    SWISSPROT       Ref Seq
  #5       P04217  NM_130786
  #7       P01023  NM_000014
  <ENS1>;
  while (<ENS1>) {
    chomp;
    my @array = split(/\t/,$_);
    my $hgnc = $array[0];
    
    if ($array[1]) {  #swissprot
      my $master = $swiss{$array[1]};
      my $dep    = $hugo{$hgnc};
      if(!defined($master) or !defined($dep)){
	$mismatch++;
      }
      else{
	XrefParser::BaseParser->add_to_xrefs($master,$hgnc,'',$hugo{$hgnc},"","",$source_id,$species_id,$count);
	$count++;
      }
      #	print "$array[1]\tSPTR\t$hgnc\tHUGO\t$hugo_id{$hgnc}\t$hugo_syn{$hgnc}\tXREF\n";
    }
    
    if ($array[2]) {
      my $master = $refseq{$array[2]};
      my $dep    = $hugo{$hgnc};
      if(!defined($master) or !defined($dep)){
	$mismatch++;
      }
      else{
	XrefParser::BaseParser->add_to_xrefs($master,$hgnc,'',$hugo{hgnc},"","",$source_id,$species_id);
	$count++;
      }
    }
  }
  close (ENS1);
  print "\t$count xrefs succesfully loaded\n";
  print "\t$mismatch xrefs failed to load\n";
}
  
      
sub new {

  my $self = {};
  bless $self, "XrefParser::HUGOParser";
  return $self;

}
 
1;
    


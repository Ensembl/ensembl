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
    print "\nUsage: HUGOParser.pm file\n\n";
    exit(1);
  }
  
  run(@ARGV);
}
  
 
sub run {
  my $self = shift if (defined(caller(1)));
  my $file = shift;
 
  my $source_id = shift;
  my $species_id = shift;

  print STDERR $file."\n";
  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
    print "source id is $source_id \n";
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
    print "species id is $species_id \n";
  }

  my $dir = dirname($file);
 
#  print STDERR $dir."\n";
  
  $xref_sth = XrefParser::BaseParser->dbi->prepare("INSERT INTO xref (accession,label,description,source_id,species_id) VALUES(?,?,?,?,?)");
  $dep_sth = XrefParser::BaseParser->dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)"); # xref1,xref2,"",hugo
    

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
    
    my $id = get_xref($hgnc, $source_id);
    if(!defined($id)){
      $xref_sth->execute($hgnc,$label,"",$source_id,$species_id);
    }
    $hugo{$hgnc} = get_xref($hgnc, $source_id);;
  }
  close ENS4;

  
  my (%swiss)  = BaseParser->get_valid_codes("uniprot",$species_id);
  my (%refseq) = BaseParser->get_valid_codes("refseq",$species_id);
  
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
	print STDERR $_."\n"; 
	print STDERR "swiss prot $array[1] -> xref $master \n";
	print STDERR "hugo number $hgnc -> xref $dep \n";
      }
      else{
	$dep_sth->execute($master, $dep,  "", $source_id);
      }
      #	print "$array[1]\tSPTR\t$hgnc\tHUGO\t$hugo_id{$hgnc}\t$hugo_syn{$hgnc}\tXREF\n";
    }
    
    if ($array[2]) {
      my $master = $refseq{$array[2]};
      my $dep    = $hugo{$hgnc};
      if(!defined($master) or !defined($dep)){
	print STDERR $_."\n"; 
	print STDERR "ref seq $array[2] -> xref $master \n";
	print STDERR "hugo number $hgnc -> xref $dep \n";
      }
      else{
#	$dep_sth->execute($master, $dep,  "", $source_id);
      }
      #	print "$array[2]\tRefSeq\t$hgnc\tHUGO\t$hugo_id{$hgnc}\t$hugo_syn{$hgnc}\tXREF\n";
    }
  }
  close (ENS1);
  
    
}


sub get_xref{
  my ($acc,$source) = @_;
  
  my $dbi = XrefParser::BaseParser->dbi;
  my $sql = "select xref_id from xref where accession = '".$acc."' and source_id = $source";
  my $sth = $dbi->prepare($sql);
    
  $sth->execute() || die $dbi->errstr;
  if(my @row = $sth->fetchrow_array()) {
    return $row[0];
  }
  return undef;
}
  
sub new {

  my $self = {};
  bless $self, "XrefParser::HUGOParser";
  return $self;

}
 
1;
    


package XrefParser::InterproParser;
  
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
    print "\nUsage: InterproParser.pm file\n\n";
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
    print "source id is $source_id \n";
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
    print "species id is $species_id \n";
  }

  my $dir = dirname($file);
 
  my %short_name;
  my %description;
  my %pfam;
    
  open (SHORT, $dir."/short_name.dat") || die "Can't open hugo interpro file $dir/short_name.dat\n";
#IPR000001       Kringle
#IPR000002       Fizzy
#IPR000003       RtnoidX_receptor

  while(<SHORT>){
    chomp;
    my ($interpro, $name) = split(/\t/,$_);
    $short_name{$interpro} = $name;
  }
  close SHORT;

  
  my $count = 0;

  open (LONG, $dir."/protein2interpro.dat") || 
    die "Can't open interpro file  $dir/protein2interpro.dat\n";
  #O00050  PF03184 166     377     IPR004875       CENP-B protein
  #O00050  PF05225 15      60      IPR007889       Helix-turn-helix, Psq
  #O00050  PF06465 448     500     IPR009463       Protein of unknown function DUF10
  #  0        1     2       3        4               5

  while (<LONG>) {
    chomp;
    my @array = split(/\t/,$_);
    my $hgnc = $array[0];
    $description{$array[4]} = $array[5];
    $pfam{$array[4]} = $array[1];
  }
  close (LONG);

  my $add_interpro_sth =  XrefParser::BaseParser->dbi->prepare
    ("INSERT INTO interpro (interpro, pfam) VALUES(?,?)");

  foreach my $interpro (keys %short_name){
    $count++;
#    print $short_name{$interpro}."\t".$interpro."\t".$description{$interpro}.
#      "\t".$pfam{$interpro}."\n";
    XrefParser::BaseParser->add_xref($interpro,'',$short_name{$interpro},
				     $description{$interpro},$source_id,$species_id);
    if(defined($pfam{$interpro})){
      $add_interpro_sth->execute($interpro,$pfam{$interpro});
    }
    else{
      print "No pfam for $interpro\n";
    }
       
  }
  print "$count xref successfully loaded.\n";
#  die "not ready yet\n";
  
  
}

sub new {

  my $self = {};
  bless $self, "XrefParser::InterproParser";
  return $self;

}
 
1;
    


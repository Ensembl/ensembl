package XrefParser::MGDParser;
 
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
    print "\nUsage: MGDParser.pm file\n\n";
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
                                                                                                                      
  my (%swiss) = XrefParser::BaseParser->get_valid_codes("uniprot",$species_id);
  my (%refseq) = XrefParser::BaseParser->get_valid_codes("refseq",$species_id);


  my $count = 0;
  my $mismatch = 0;

  open(FILE,"<". $file) || die "could not open file $file";
  while(my $line = <FILE>){
    chomp $line;
    my ($key,$label,$sps) = (split("\t",$line))[0,1,6];
    my @sp = split(/\s/,$sps); 
    foreach my $value (@sp){
      if(defined($value) and $value and defined($swiss{$value})){
	XrefParser::BaseParser->add_to_xrefs($swiss{$value},$key,'',$label,"",$source_id,$species_id);
	$count++;
      }
      elsif(defined($value) and $value and defined($refseq{$value})){
#	print STDERR "ARSE\t".$key."\t".$value."\t".$label."\n";	
	$mismatch++;
      }
    }
  }
  close FILE;
  print "$count xrefs succesfully loaded\n";
  print "$mismatch xrefs failed to load\n";
     



}                                                                                                                     

sub new {

  my $self = {};
  bless $self, "XrefParser::MGDParser";
  return $self;

}
 
1;
    

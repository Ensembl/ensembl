package XrefParser::MGDParser;
 
use strict;
use POSIX qw(strftime);
use File::Basename;
 
use BaseParser;
 
use vars qw(@ISA);
@ISA = qw(BaseParser);

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

  $xref_sth = BaseParser->dbi->prepare("INSERT INTO xref (accession,label,description,source_id,species_id) VALUES(?,?,?,?,?)");
                                                                                                                      
  $dep_sth = BaseParser->dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");

#  my $mrk_swiss = "MRK_SwissProt_TrEMBL.rpt";
#  my $mrk_ll    = "MRK_LocusLink.rpt";

#  my $dir =  dirname($file);  
  
  if(!defined($source_id)){
    $source_id = BaseParser->get_source_id_for_filename($file);
    print "source id is $source_id \n";
  }
  if(!defined($species_id)){
    $species_id = BaseParser->get_species_id_for_filename($file);
    print "species id is $species_id \n";
  }                                                                                                                      
                                                                                                                      
  my (%swiss) = BaseParser->get_valid_codes("uniprot",$species_id);
  my (%refseq) = BaseParser->get_valid_codes("refseq",$species_id);


  print %swiss."\n";


  open(FILE,"<". $file) || die "could not open file $file";
  while(my $line = <FILE>){
    chomp $line;
    my ($key,$label,$sps) = (split("\t",$line))[0,1,6];
    my @sp = split(/\s/,$sps); 
    foreach my $value (@sp){
#      print $value."\n";
      if(defined($value) and $value and defined($swiss{$value})){
	add_to_xrefs($swiss{$value},$key,$label,$source_id,$species_id);

#	print $key."\t".$value."\t".$label."\n";
      }
      elsif(defined($value) and $value and defined($refseq{$value})){
	print STDERR "ARSE\t".$key."\t".$value."\t".$label."\n";	
      }
    }
  }
  close FILE;
     



}                                                                                                                     

sub get_xref{
  my ($acc,$source) = @_;
 
  my $dbi =BaseParser->dbi;
  my $sql = "select xref_id from xref where accession = '".$acc."' and source_id = $source";
  my $sth = $dbi->prepare($sql);
   
  $sth->execute() || die $dbi->errstr;
  if(my @row = $sth->fetchrow_array()) {
    return $row[0];
  }
  return undef;
}
 
sub add_to_xrefs{
  my ($master_xref,$mgi,$label,$source_id,$species_id) = @_;
 
  my $dependent_id = get_xref($mgi, $source_id);
  if(!defined($dependent_id)){
    $xref_sth->execute($mgi,$label,"",$source_id,$species_id);
  }
  $dependent_id = get_xref($mgi, $source_id);
  $dep_sth->execute($master_xref, $dependent_id,  "", $source_id);
 
}

package GoParser;

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
    print "\nUsage: GoParser.pm file\n\n";
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

  open(GO,"<".$file) || die "Could not open $file\n";

  while (<GO>) {
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
	add_to_xrefs($refseq{$array[1]},$array[4],$array[6],$source_id,$species_id);
#	print "$array[1]\tSPTR\t$array[4]\tGO\t$array[6]\t$array[9]\tXREF\n";
      }
    }
    elsif($array[0] =~ /UniProt/){
      if($swiss{$array[1]}){
	add_to_xrefs($swiss{$array[1]},$array[4],$array[6],$source_id,$species_id);
#	print "$array[1]\tSPTR\t$array[4]\tGO\t$array[6]\t$array[9]\tXREF\n";
      }
    }
    else{
      print STDERR "unknown type ".$array[0]."\n";
    }
  }
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
  my ($master_xref,$go,$linkage,$source_id,$species_id) = @_;

  my $dependent_id = get_xref($go, $source_id);
  if(!defined($dependent_id)){
    $xref_sth->execute($go,$go,"",$source_id,$species_id);
  }
  $dependent_id = get_xref($go, $source_id);
  $dep_sth->execute($master_xref, $dependent_id,  $linkage, $source_id);

}

package XrefParser::WormPepParser;

use strict;
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
    print "\nUsage: WormPepParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run(@ARGV);
}

sub run {

  my $self = shift if (defined(caller(1)));
  my $file = shift;

  my $source_id = shift;
  my $species_id = shift;

  print STDERR "WORMPep source = $source_id\tspecies = $species_id\n";
  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }
  my %worm = %{XrefParser::BaseParser->get_valid_xrefs_for_dependencies
	      ('wormbase_transcript','Uniprot/SPTREMBL','RefSeq_peptide',
	       'Uniprot/SWISSPROT')}; 
  my $sql = "update xref set accession =\'?\' where xref_id=?";
  my $dbi =  XrefParser::BaseParser->dbi();
  my $sth = $dbi->prepare($sql);


  my $count  = 0;

  open(PEP,"<".$file) || die "Could not open $file\n";

  while (<PEP>) {
    my ($transcript, $wb, $swiss)  = (split(/\t/,substr($_,1)))[0,1,5];
    
    if($swiss =~ /SW:(.*)/){
      $swiss = $1;
    }
#    else{
#      print "ERR:".$_;
#    }
    if(length($wb) < 3){
      print "ERRR:".$_;
    }
    
    #Is the transcript different from the gene
    my $diff =0;
    my $gene;
    if($transcript =~ /(\S+\.\d+)/){
      $gene = $1;
      if($gene ne $transcript){
	$diff=1;
      }
    }
    else{
      die "Gene format not recognised $transcript\n";
    }

    my $exists =0;
    # if gene stored as transcript so change this
    if(defined($worm{$gene}) and !defined($worm{$transcript})){
      # change accesion to transcript name instead of gene
      $sth->execute($transcript, $worm{$gene}) || die $dbi->errstr;
      print "changing $gene (".$worm{$gene}.") to $transcript\n";
    }
    # if no record exists for this 
    elsif(!defined($worm{$gene}) and !defined($worm{$transcript})){
      print "Could not find: $transcript  ????";
      if($swiss){
	print " on sw: $swiss";
      }
      print "\n";
    }
      

#    print $gene,"\t",$swiss,"\t",$wb,"\n";
  }
#  die "Just testing at the moment\n";

}

  
sub new {

  my $self = {};
  bless $self, "XrefParser::WormPepParser";
  return $self;

}

1;


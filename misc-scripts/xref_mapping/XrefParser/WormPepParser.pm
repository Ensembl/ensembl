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

  my $worm_source_id = XrefParser::BaseParser->get_source_id_for_source_name('wormbase_transcript');

  my (%worm)  =  %{XrefParser::BaseParser->get_valid_codes("wormbase_transcript",$species_id)};

  my (%swiss)  =  %{XrefParser::BaseParser->get_valid_codes("Uniprot",$species_id)};

  my $sql = "update xref set accession =? where xref_id=?";
  my $dbi =  XrefParser::BaseParser->dbi();
  my $sth = $dbi->prepare($sql);


  my $sql2 = "select x2.accession, x2.xref_id ";
  $sql2   .= "from dependent_xref d, xref x1, xref x2 ";
  $sql2   .= "where d.master_xref_id = x1.xref_id and ";
  $sql2   .= "      d.dependent_xref_id = x2.xref_id and ";
  $sql2   .= "      x2.source_id = $worm_source_id and ";
  $sql2   .= "      x1.xref_id = ? and ";
  $sql2   .= "      x2.accession = ?";
  my $sth2 = $dbi->prepare($sql2);


  my $sql3 = 'delete from dependent_xref where dependent_xref.master_xref_id=? and dependent_xref.dependent_xref_id=?'; 
  my $sth3 = $dbi->prepare($sql3);

  open(PEP,"<".$file) || die "Could not open $file\n";

  while (<PEP>) {
    my ($transcript, $wb, $swiss_ref)  = (split(/\t/,substr($_,1)))[0,1,5];
    my $swiss_xref;
    if($swiss_ref =~ /SW:(.*)/){
      $swiss_ref = $1;
      if(defined($swiss{$swiss_ref})){
	$swiss_xref = $swiss{$swiss_ref};
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

	$sth2->execute($swiss_xref, $gene) || die $dbi->errstr;
	(my $gene_acc, my $gene_xref) = $sth2->fetchrow_array();

	$sth2->execute($swiss_xref, $transcript) || die $dbi->errstr;
	(my $tran_acc, my $tran_xref) =  $sth2->fetchrow_array();

	my $create = 1;
	if(defined($tran_xref)){ #okay
	  $create = 0;
	}
	elsif(defined($gene_xref)){
	  #need to delete dependency
	  #then add new one with correct name
	  $sth3->execute($swiss_xref, $gene_xref) || die $dbi->errstr;
	  print "removing $swiss_ref -> $gene : ";
	}
	if($create){
	  XrefParser::BaseParser->add_to_xrefs($swiss_xref,$transcript,'',$transcript,"","",$worm_source_id,$species_id);	  
	  print "adding $swiss_ref -> $transcript\n";
	}
	
      }
      
    }
  }
  
}

  
sub new {

  my $self = {};
  bless $self, "XrefParser::WormPepParser";
  return $self;

}

1;


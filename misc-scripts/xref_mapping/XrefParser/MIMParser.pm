package XrefParser::MIMParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);



# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: MIMParser.pm file <source_id> <species_id>\n\n";
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


  my (%genename2xref) = gene_name_2_xref_from_hugo();

  my $count =0;
  my $mismatch=0;

  open(MIM,"<".$file) || die "Could not open $file\n";

  while (<MIM>) {
    chomp;
    my ($desc, $gene_names, $mim_number, $loc) = split (/\|/,$_);

    my $xref =0;
    foreach my $gene (split(/\, */,$gene_names)){
      if(defined($genename2xref{$gene})){
	$xref = $genename2xref{$gene};
      }
    }
    if($xref){
       XrefParser::BaseParser->add_to_xrefs($xref,$mim_number,'',$mim_number,$desc,'',$source_id,$species_id);  
      $count++;
    }
    else{
       $mismatch++;
    }
  }
  print "\t$count succesfull xrefs loaded\n";
  print "\t$mismatch FAILED xrefs\n";
}

sub gene_name_2_xref_from_hugo{
  my %gene_name2xref;
  
  my $dbi = XrefParser::BaseParser->dbi();


  my $source_id_for_hugo=0;
  my $sth2 = $dbi->prepare("select * from source where name like ?");
  $sth2->execute('HUGO') || die $dbi->errstr;
  while(my @row = $sth2->fetchrow_array()) {
    $source_id_for_hugo = $row[0];
  }
  if(! $source_id_for_hugo){
    die "Could not find source id for HUGO.\n";
  }


  my @source_list=();

  $sth2->execute('Refseq') || die $dbi->errstr;
  while(my @row = $sth2->fetchrow_array()) {
    push  @source_list, $row[0];
  }

  $sth2->execute('Uniprot') || die $dbi->errstr;
  while(my @row = $sth2->fetchrow_array()) {
    push  @source_list, $row[0];
  }


  $sth2->finish;


  my $sql = "select y.label, x.xref_id ";
  $sql   .= "  from dependent_xref d, xref x, xref y ";
  $sql   .= "  where d.source_id = $source_id_for_hugo and ";
  $sql   .= "        x.xref_id = d.master_xref_id and ";
  $sql   .= "        y.xref_id = d.dependent_xref_id and ";
  $sql   .= "        x.source_id = ?";
  
    
  my $sth = $dbi->prepare($sql);

  foreach my $id (@source_list){
    $sth->execute($id) || die $dbi->errstr;
    while(my @row = $sth->fetchrow_array()) {
      my $gene_name = $row[0];
      my $xref = $row[1];
      $gene_name2xref{$gene_name} = $xref;
    }
  } 
  $sth->finish;
  return %gene_name2xref;
}
  


sub new {

  my $self = {};
  bless $self, "XrefParser::MIMParser";
  return $self;

}
 
1;

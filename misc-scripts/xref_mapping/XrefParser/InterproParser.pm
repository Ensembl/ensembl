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
    print "\nUsage: InterproParser.pm file <source_id> <species_id>\n\n";
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

  my $add_interpro_sth =  XrefParser::BaseParser->dbi->prepare
    ("INSERT INTO interpro (interpro, pfam) VALUES(?,?)");

  my $get_interpro_sth =  XrefParser::BaseParser->dbi->prepare
    ("SELECT interpro from interpro where interpro = ? and pfam = ?");
 
  my $add_xref_sth = XrefParser::BaseParser->dbi->prepare
    ("INSERT INTO xref (accession,version,label,description,source_id,species_id) VALUES(?,?,?,?,?,?)");
  
  my $get_xref_sth = XrefParser::BaseParser->dbi->prepare
    ("select xref_id from xref where accession = ? and source_id = ?");


  my $dir = dirname($file);
                                                                                                                         
  my %short_name;
  my %description;
  my %pfam;
     
  open (XML, $dir."/interpro.xml") || die "Can't open hugo interpro file $dir/interpro.xml\n";
  #<interpro id="IPR001023" type="Family" short_name="Hsp70" protein_count="1556">
  #    <name>Heat shock protein Hsp70</name>
  #     <db_xref protein_count="18" db="PFAM" dbkey="PF01278" name="Omptin" />
  #      <db_xref protein_count="344" db="TIGRFAMs" dbkey="TIGR00099" name="Cof-subfamily" />
  
  my %count;
  local $/ = "</interpro>";
 

  my $last = "";
  my $i =0;
  while (<XML>) {

    my $interpro;
    my $short_name;
    my $name;
    
    ($interpro) = $_ =~ /interpro id\=\"(\S+)\"/;
    ($short_name) = $_ =~ /short_name\=\"(\S+)\"/;
    ($name) = $_ =~ /\<name\>(.*)\<\/name\>/;
	  
    if($interpro){
#      print $interpro."\n";
      if(!get_xref($get_xref_sth, $interpro, $source_id)){
	$count{INTERPRO}++;
	$add_xref_sth->execute($interpro,'',$short_name, $name,$source_id,$species_id)
	  || die "Problem adding ".$interpro."\n";
      }

      while( /db="(PROSITE|PFAM|PRINTS|PREFILE|PROFILE|TIGRFAMs)"\s+dbkey="(\S+)"/cgm ) {
	my ( $db_type, $id ) =  ( $1, $2 );
#	print $db_type."\t".$id."\n";
	if(!get_xref($get_interpro_sth, $interpro,$id)){
	  $add_interpro_sth->execute($interpro,$id);
	  $count{$db_type}++;
	}
      }
    }
  }

  close (LONG);
#  die "\n";
  for my $db ( keys %count ) {
    print "\t".$count{$db}." $db loaded.\n";
  }
  
}

sub get_xref{
  my ($get_xref_sth, $acc, $source) = @_;

  $get_xref_sth->execute($acc, $source) || die "FAILED $acc  $source\n";
  if(my @row = $get_xref_sth->fetchrow_array()) {
    return $row[0];
  }   
  return 0;
}

sub new {

  my $self = {};
  bless $self, "XrefParser::InterproParser";
  return $self;

}
 
1;
    


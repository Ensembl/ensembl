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
  
  my $count  = 0;
  my $count2 = 0;
  my $count3 = 0;
  local $/ = "</interpro>";
 

  my $last = "";
  my $i =0;
  while (<XML>) {

    my $interpro;
    my $short_name;
    my $pfam;
    my $tigr;
    my $name;
    
    ($interpro) = $_ =~ /interpro id\=\"(\S+)\"/;
    ($short_name) = $_ =~ /short_name\=\"(\S+)\"/;
    ($name) = $_ =~ /\<name\>(.*)\<\/name\>/;
    ($pfam) = $_ =~ /db\=\"PFAM\".*dbkey\=\"(\S+)\"/;
    ($tigr) = $_ =~ /db\=\"TIGRFAMs\".*dbkey\=\"(\S+)\"/;
    

#    print "#########################################################\n$interpro\n$name\n$short_name\n";
#    $i++;
#    if($i > 10){
#      die "first ten done";
#    }
    if($interpro){
      if(!get_xref($get_xref_sth, $interpro, $source_id)){
	$count++;
	$add_xref_sth->execute($interpro,'',$short_name, $name,$source_id,$species_id)
	  || die "Problem adding ".$interpro."\n";
      }
      if($pfam){
	#      print "PFAM $pfam\n";
	if(!get_xref($get_interpro_sth, $interpro,$pfam)){
	  $add_interpro_sth->execute($interpro,$pfam);
	  $count2++;
	}
      }  
      if($tigr){
	#     print "TIGR $tigr\n";
	if(!get_xref($get_interpro_sth, $interpro,$tigr)){
	  $add_interpro_sth->execute($interpro,$tigr);
	  $count3++;
	}
      }  
    }
  }
  close (LONG);
       
  print "$count xref successfully loaded.\n";
  print "$count2 interpro/pfam relationships added\n";
  print "$count3 interpro/tigr relationships added\n";
#  die "not ready yet\n";
  
  
}

sub get_xref{
  my ($get_xref_sth, $acc, $source) = @_;

  $get_xref_sth->execute($acc, $source) || die "FAILED $acc  $source\n";
  if(my @row = $get_xref_sth->fetchrow_array()) {
#    print "FOUND $acc\n";
    return $row[0];
  }   
#  print "UNKOWN $acc";
  return 0;
}

sub new {

  my $self = {};
  bless $self, "XrefParser::InterproParser";
  return $self;

}
 
1;
    


package XrefParser::ncRNAParser;
 
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
    print "\nUsage: ncRNAParser.pm file <source_id> <species_id>\n\n";
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

  my %name_2_source_id=();
  my $added=0;

  if(!open(FILE,"<". $file)){
    print  "ERROR: Could not open file $file\n";
    return 1;
  }
  while(my $line = <FILE>){
    chomp $line;
    my ($gene_id,$transcript_id,$source_name,$acc,$display_label,$full_description)
      = split("\t",$line);

    #trim the description.
    my ($description,$junk) = split("[[]Source:",$full_description);
    if(!defined($name_2_source_id{$source_name})){
      my $tmp = $self->get_source_id_for_source_name($source_name);
      if(!$tmp){
	die("Could not get source_id for $source_name\n");
      }
      $name_2_source_id{$source_name} = $tmp;
    }
    my $xref_id = $self->get_xref($acc,$name_2_source_id{$source_name});
    if(!defined($xref_id)){
      $xref_id = $self->add_xref($acc,"",$display_label,$description,$name_2_source_id{$source_name}, $species_id);
      $added++;
    }
    $self->add_direct_xref($xref_id, $transcript_id, "transcript", "") if (defined($transcript_id));    
    $self->add_direct_xref($xref_id, $gene_id, "gene", "")             if (defined($gene_id)); 
  }
  close FILE;

  print "Added $added Xrefs for ncRNAs\n";
  return 0;
}

sub new {

  my $self = {};
  bless $self, "XrefParser::ncRNAParser";
  return $self;

}


1;

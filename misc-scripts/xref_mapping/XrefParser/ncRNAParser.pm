package XrefParser::ncRNAParser;
 
use strict;
use POSIX qw(strftime);
use File::Basename;
 
use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: ncRNAParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0]);

}


sub run {
  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  print "source_id = $source_id, species= $species_id, file = $file\n";


  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }                                                

  my %name_2_source_id=();
  my $added=0;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open file $file\n";
    return 1;
  }

  while ( my $line = $file_io->getline() ) {
    chomp $line;
    my ($gene_id,$transcript_id,$source_name,$acc,$display_label,$full_description, $status)
      = split("\t",$line);

    #trim the description.
    my ($description,$junk) = split("[[]Source:",$full_description);
    if($source_name eq "miRNA_Registry"){
      if($status eq "KNOWN"){
	$source_name = "miRBase";
      }
      else{
	$source_name = "miRBase_predicted";
      }
    }
    if(!defined($name_2_source_id{$source_name})){
      my $tmp = $self->get_source_id_for_source_name($source_name);
      if(!$tmp){
	die("Could not get source_id for $source_name\n");
      }
      $name_2_source_id{$source_name} = $tmp;
    }
    my $xref_id = $self->get_xref($acc,$name_2_source_id{$source_name}, $species_id);
    if(!defined($xref_id)){
      $xref_id = $self->add_xref($acc,"",$display_label,$description,$name_2_source_id{$source_name}, $species_id,"DIRECT");
      $added++;
    }
    $self->add_direct_xref($xref_id, $transcript_id, "Transcript", "") if (defined($transcript_id));    

    #just add to the transcript ONLY as the check at the end will move all
    #the those mapped to the transcript to the genes anyway due to the
    #biomart check
    #    $self->add_direct_xref($xref_id, $gene_id, "Gene", "")             if (defined($gene_id)); 
  }

  $file_io->close();

  print "Added $added Xrefs for ncRNAs\n" if($verbose);
  return 0;
}

1;

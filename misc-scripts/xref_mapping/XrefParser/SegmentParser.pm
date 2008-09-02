package XrefParser::SegmentParser;
 
use strict;
use POSIX qw(strftime);
use File::Basename;
 
use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: SegmentParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0]);

}

sub run {
  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $files_ref  = shift;
  my $rel_file   = shift;
  my $verbose = shift;

  my $file = @{$files_ref}[0];


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
    my ($gene_id,$transcript_id,$source_name,$acc,$display_label,$description, $status)
      = split("\t",$line);

      $source_name =~ tr/-/_/;  # No minuses allowed, change to underscore.

    if(!defined($name_2_source_id{$source_name})){
      my $tmp = $self->get_source_id_for_source_name($source_name);
      if(!$tmp){
	die("Could not get source_id for $source_name\n");
      }
      $name_2_source_id{$source_name} = $tmp;
    }
    my $xref_id = $self->get_xref($acc,$name_2_source_id{$source_name}, $species_id);
    if(!defined($xref_id)){
      $xref_id = $self->add_xref($acc,"",$display_label,$description,$name_2_source_id{$source_name}, $species_id);
      $added++;
    }
    $self->add_direct_xref($xref_id, $transcript_id, "Transcript", "") if (defined($transcript_id));    

    #just add to the transcript ONLY as the check at the end will move all
    #the those mapped to the transcript to the genes anyway due to the
    #biomart check
  }

  $file_io->close();

  print "Added $added Xrefs for Gene segments\n" if($verbose);
  return 0;
}

1;

package XrefParser::SegmentParser;
 
use strict;
use POSIX qw(strftime);
use File::Basename;
 
use base qw( XrefParser::BaseParser );

sub run {
  my $self = shift;
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

  my $added=0;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open file $file\n";
    return 1;
  }

  while ( my $line = $file_io->getline() ) {
    chomp $line;
    my ($gene_id,$acc,$description)  = split("\t",$line);

    my $xref_id = $self->get_xref($acc,$source_id, $species_id);
    if(!defined($xref_id)){
      $xref_id = $self->add_xref($acc,"",$acc,$description,$source_id, $species_id, "DIRECT");
      $added++;
    }
#    print "$acc, $xref_id, $gene_id\n";
    $self->add_direct_xref($xref_id, $gene_id, "Gene", "");

  }

  $file_io->close();

  print "Added $added Xrefs for Gene segments\n" if($verbose);
  return 0;
}

1;

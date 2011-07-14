package XrefParser::VBRNADescriptionParser;

use strict;
use POSIX qw(strftime);
use File::Basename;
use base qw( XrefParser::BaseParser );

# Parse the RNA description file
#
# AGAP013564      5.8S ribosomal RNA [Source: RFAM 9.0]
# AGAP013555      Small nucleolar RNA U3 [Source: RFAM 9.0]
# AGAP013556      U6 spliceosomal RNA [Source: RFAM 9.0]
# etc.

sub run {

  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  print "source_id = $source_id, species= $species_id, file = $file\n" if($verbose);

  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }

  my $added = 0;
  my $count = 0;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open file $file\n";
    return 1;
  }

  while ( my $line = $file_io->getline() ) {
    chomp $line;
    my ($gene_id, $description) = split("\t",$line);  #and use the gene_id as accession

    my $xref_id = $self->get_xref($gene_id,$source_id, $species_id);
    if(!defined($xref_id)){
      $xref_id = $self->add_xref($gene_id,"", $gene_id, $description, $source_id, $species_id, "DIRECT");
      $count++;
    }
    if(defined($gene_id) and $gene_id ne "-"){
      $self->add_direct_xref($xref_id, $gene_id, "Gene", "") ;
      $added++;
    }	
  }

  $file_io->close();

  print "Added $count xrefs and $added Direct xrefs to genes for VBRNADescription\n" if($verbose);
  return 0;
}

1;

  

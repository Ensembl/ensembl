package XrefParser::AnoXcelParser;

use strict;
use POSIX qw(strftime);
use File::Basename;
use base qw( XrefParser::BaseParser );

# Parse the external description file
#
# Protein         mRNA            Location                Gene stable_id
# AGAP004677-PB   AGAP004677-RB   2L:157496-159356:-1     AGAP004677
# AGAP004677-PA   AGAP004677-RA   2L:157496-181213:-1     AGAP004677
# AGAP004678-PA   AGAP004678-RA   2L:203866-204956:1      AGAP004678
#...

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
    my ($protein_id, $mRNA_id, $loc, $gene_id) = split("\t",$line);  #and use the gene_id as accession

    my $xref_id = $self->get_xref($gene_id,$source_id, $species_id);
    if(!defined($xref_id)){
      $xref_id = $self->add_xref($gene_id,"", $gene_id, $gene_id, $source_id, $species_id, "DIRECT");
      $count++;
    }
    if(defined($gene_id) and $gene_id ne "-"){
      $self->add_direct_xref($xref_id, $gene_id, "Gene", "") ;
      $added++;
    }
  }

  $file_io->close();

  print "Added $count xrefs and $added Direct xrefs to genes for AnoXcel\n" if($verbose);
  return 0;
}
1;


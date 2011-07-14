package XrefParser::ImmunoDBParser;

use strict;
use POSIX qw(strftime);
use File::Basename; 
use base qw( XrefParser::BaseParser );

sub run {
  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  print "source_id = $source_id, species = $species_id, file = $file\n" if($verbose);

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

    my ($SPECIES,$gene_id, $acc, $family, $subfamily, $description) = split(",",$line);

    my $full_description = $description."($family)" ;
    if ($subfamily ne $family) { $full_description .= ", subfamily $subfamily" ;}
    #$subfamily ~= /1-3-beta-D/1,3-beta-D/ ;

    my $xref_id = $self->get_xref($acc,$source_id, $species_id);

    if(!defined($xref_id)){
      $xref_id = $self->add_xref($acc,"", $acc, $full_description, $source_id, $species_id, "DIRECT");
      $count++;
    }
    if(defined($gene_id) and $gene_id ne "-"){
      $self->add_direct_xref($xref_id, $gene_id, "Gene", "") ;
      $added++;
    }	
  }

  $file_io->close();

  print "Added $count xrefs and $added Direct xrefs to genes for ImmunoDB\n" if($verbose);
  return 0;
}
1;

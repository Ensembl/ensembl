package XrefParser::VBRNADescriptionParser;

use strict;
use warnings;
use Carp;
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

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  print "source_id = $source_id, species= $species_id, file = $file\n" if($verbose);

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

  

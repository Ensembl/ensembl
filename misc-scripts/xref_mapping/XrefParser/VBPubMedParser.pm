package XrefParser::VBPubMedParser;

use strict;
use POSIX qw(strftime);
use File::Basename;
use base qw( XrefParser::BaseParser );

# Parse the external description file
#
#PubMed ID      Stable ID       Feature Gene_ID (?)     Origin
#1354853	AAEL009742	gene	AAEL009742	Inferred from UniProt entry ABDA_AEDAE (P29552)
#1961751	AAEL006563	gene	AAEL006563	Inferred from UniProt entry VCP_AEDAE (P42660)
#2052024	AAEL006424	gene	AAEL006424	Inferred from UniProt entry ALL2_AEDAE (P18153)


if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: VBPubMed.pm file <source_id> <species_id>\n\n";
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
    if ($line !~ /^#/) {

      chomp $line;
      my ($PubMed_id, $gene_id, $rien, $name, $origin) = split("\t",$line);  #and use the gene_id as accession
      my $descr_full = "PubMed ID $PubMed_id - $origin\n" ;

      print "PMID: $PubMed_id, GID: $gene_id, NONE: $rien, NOM: $name, ORI: $origin\n" ;

      my $xref_id = $self->get_xref($gene_id,$source_id, $species_id);
      if(!defined($xref_id)){
	$xref_id = $self->add_xref($gene_id,"", $gene_id, $descr_full, $source_id, $species_id, "DIRECT");
	$count++;
      }
      if(defined($gene_id) and $gene_id ne "-"){
	$self->add_direct_xref($xref_id, $gene_id, "Gene", "") ;
	$added++;
      }	
    }

    $file_io->close();

    print "Added $count xrefs and $added Direct xrefs to genes for VBPubMed\n" if($verbose);
    return 0;
  }
}

1;



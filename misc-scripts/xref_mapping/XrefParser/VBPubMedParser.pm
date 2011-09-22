package XrefParser::VBPubMedParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;
use base qw( XrefParser::BaseParser );

# Parse the external description file
#
#PubMed ID      Stable ID       Feature Gene_ID (?)     Origin
#1354853	AAEL009742	gene	AAEL009742	Inferred from UniProt entry ABDA_AEDAE (P29552)
#1961751	AAEL006563	gene	AAEL006563	Inferred from UniProt entry VCP_AEDAE (P42660)
#2052024	AAEL006424	gene	AAEL006424	Inferred from UniProt entry ALL2_AEDAE (P18153)

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
    if ($line !~ /^#/) {

      chomp $line;
      my ($PubMed_id, $gene_id, $rien, $name, $origin) = split("\t",$line);  #and use the gene_id as accession
      my $descr_full = "PubMed ID $PubMed_id - $origin" ;
      #print STDERR "PMID: $PubMed_id, GID: $gene_id, NONE: $rien, NOM: $name, ORI: $origin\n" ;

      my $xref_id = $self->get_xref($gene_id,$source_id, $species_id);
      if(!defined($xref_id)){
	$xref_id = $self->add_xref({ acc        => $gene_id,
				     label      => $gene_id,
				     desc       => $descr_full,
				     source_id  => $source_id,
				     species_id => $species_id,
				     info_type  => "DIRECT"} );
	$count++;
      }
      if(defined($gene_id) and $gene_id ne "-"){
	$self->add_direct_xref($xref_id, $gene_id, "Gene", "") ;
	$added++;
      }
    }
  }
  $file_io->close();

  print "Added $count xrefs and $added Direct xrefs to genes for VBPubMed\n" if($verbose);
  return 0;
}

1;



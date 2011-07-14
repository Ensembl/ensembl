package XrefParser::SGDParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

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
  
  my $gene_source_id = XrefParser::BaseParser->get_source_id_for_source_name("SGD_GENE");
  my $transcript_source_id = XrefParser::BaseParser->get_source_id_for_source_name("SGD_TRANSCRIPT");

  my $sgd_io = $self->get_filehandle($file);

  if ( !defined $sgd_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $xref_count =0;
  my $syn_count =0;

  while ( $_ = $sgd_io->getline() ) {

    chomp;

    if ($_ =~ /^([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t[^\t]+\t[^\t]*\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]*\t[^\t]+\t[^\t]+\t([^\t]*)$/) {
	my ($sgd_id, $biotype, $status, $orf_name, $locus_name, $alias_name, $desc) = ($1,$2,$3,$4,$5,$6,$7);
	
	# parse the lines corresponding to the gene entries
	# and filter out lines corresponding to the CDS for example
	
	if ($biotype =~ /ORF|.+RNA|transposable_element_gene|pseudogene/) {

	    if ($verbose) {	    
		print STDERR "parsing line for biotype, $biotype\n";
	    
		print STDERR "sgd_id, biotype, status, orf_name, locus_name, alias_name, $sgd_id, $biotype, $status, $orf_name, $locus_name, $alias_name\n";
		print STDERR "desc: $desc\n";
	    }

	    if (!defined $locus_name || ($locus_name eq "")) {
		
		if (!defined $orf_name || ($orf_name eq "")) {
		    print STDERR "can't assign the orf_name as the locus_name!\n";
		}
		else {
		    if ($verbose) {
			print STDERR "assigning the orf_name as the locus_name(ie the gene_name)\n";
		    }
		    $locus_name = $orf_name;
		}
	    }
	    
	    my (@syn) = split(/\|/,$alias_name);

	    my $gene_xref_id = $self->add_xref($sgd_id,"",$locus_name,$desc,$gene_source_id,$species_id,"DIRECT");
	    $self->add_direct_xref($gene_xref_id, $orf_name, "Gene", "DIRECT");

	    my $transcript_xref_id = $self->add_xref($sgd_id,"",$locus_name,$desc,$transcript_source_id,$species_id,"DIRECT");
	    $self->add_direct_xref($transcript_xref_id, $orf_name, "Transcript", "DIRECT");

	    $xref_count++;

	    foreach my $synonym (@syn){
		if ($verbose) {
		    print STDERR "adding synonym, $synonym\n";
		}
		$self->add_to_syn($sgd_id, $gene_source_id, $synonym, $species_id);
		$syn_count++;
	    }
	    
	}
	else {
	    if ($verbose) {
		print STDERR "filtering biotype, $biotype\n";
	    }
	}
    }
    else {
	if ($verbose) {
	    print STDERR "failed to parse line, $_\n\n";
	}
    }
  }

  $sgd_io->close();

  print $xref_count." SGD Xrefs added with $syn_count synonyms\n" if($verbose);
  return 0; #successful
}

1;

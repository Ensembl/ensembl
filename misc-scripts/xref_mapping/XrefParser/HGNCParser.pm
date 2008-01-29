package XrefParser::HGNCParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

my $xref_sth ;
my $dep_sth;
my $syn_sth;

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: HGNCParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run(@ARGV);
}

sub run {

  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $file = shift;

  print STDERR "source = $source_id\tspecies = $species_id\n";
  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }


  my $hgnc_refseq_manual = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","refseq_manual");
  if(!defined($hgnc_refseq_manual)){
    die  "Could not get source id for HGNC with priority description of refseq_manual\n";
  }
  my $hgnc_refseq_mapped = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","refseq_mapped");
  if(!defined($hgnc_refseq_mapped)){
    die  "Could not get source id for HGNC with priority description of refseq_mapped\n";
  }

  my $hgnc_entrezgene_manual  = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","entrezgene_manual");
  if(!defined($hgnc_entrezgene_manual)){
    die  "Could not get source id for HGNC with priority description of entrezgene_manual\n";
  }
  my $hgnc_entrezgene_mapped  = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","entrezgene_mapped");
  if(!defined($hgnc_entrezgene_mapped)){
    die  "Could not get source id for HGNC with priority description of entrezgene_mapped\n";
  }

#  my (%swiss)  =  %{XrefParser::BaseParser->get_valid_codes("uniprot",$species_id)};
  my (%refseq) =  %{XrefParser::BaseParser->get_valid_codes("refseq",$species_id)};
  my @list;
  push @list, "refseq_peptide";
  push @list, "refseq_dna";
  my (%entrezgene) = %{XrefParser::BaseParser->get_valid_xrefs_for_dependencies("EntrezGene",@list)};

  my $swiss_count = 0;
  my $refseq_count = 0;
  my $entrezgene_count = 0;
  my $mismatch = 0;

  my $hugo_io = $self->get_filehandle($file);

  if ( !defined $hugo_io ) {
    print "ERROR: Can't open HGNC file $file\n";
    return 1;
  }

  $_ = $hugo_io->getline();

  while ( $_ = $hugo_io->getline() ) {

    chomp;

    # 0 HGNC ID	           # primary accession
    # 1 Approved Symbol    # label
    # 2 Approved Name      # description
    # 3 Previous Symbols   # synonyms
    # 4 Aliases            # aliases
    # 5 entrezgene ID   manually curated
    # 6 RefSeq ID       manually curated
    # 7 entrezgene ID   mapped
    # 8 RefSeq ID       mapped

    my @array = split(/\t/,$_);

    # Use the RefSeq if available as this is manually curated
    # If no RefSeq, use the Swissprot instead

    my $seen = 0;
    if ($array[6]) {             # RefSeq
      if(defined($refseq{$array[6]})){
	$seen = 1;
	$refseq_count++;
	XrefParser::BaseParser->add_to_xrefs($refseq{$array[6]}, $array[0], '', $array[1], $array[2], "", $hgnc_refseq_manual, $species_id);

	if (defined($array[3])) {     # dead name, add to synonym
	  my @array2 = split(',\s*', $array[3]);
	  foreach my $arr (@array2){
	    XrefParser::BaseParser->add_to_syn($array[0], $hgnc_refseq_manual, $arr);
	  }
	}

	if (defined($array[4])) {     # alias, add to synonym
	  my @array2 = split(',\s*', $array[4]);
	  foreach my $arr (@array2){
	    XrefParser::BaseParser->add_to_syn($array[0], $hgnc_refseq_manual, $arr);
	  }
	}
      }
    }
    if ($array[8]) {             # RefSeq
      if(defined($refseq{$array[8]})){
	$seen = 1;
	$refseq_count++;
	XrefParser::BaseParser->add_to_xrefs($refseq{$array[8]}, $array[0], '', $array[1], $array[2], "", $hgnc_refseq_mapped, $species_id);

	if (defined($array[3])) {     # dead name, add to synonym
	  my @array2 = split(',\s*', $array[3]);
	  foreach my $arr (@array2){
	    XrefParser::BaseParser->add_to_syn($array[0], $hgnc_refseq_mapped, $arr);
	  }
	}

	if (defined($array[4])) {     # alias, add to synonym
	  my @array2 = split(',\s*', $array[4]);
	  foreach my $arr (@array2){
	    XrefParser::BaseParser->add_to_syn($array[0], $hgnc_refseq_mapped, $arr);
	  }
	}
      }
    }

    if(defined($array[5])){
      if(defined($entrezgene{$array[5]})){
	$seen = 1;
	XrefParser::BaseParser->add_to_xrefs($entrezgene{$array[5]}, $array[0], '', 
					     $array[1], $array[2], "", $hgnc_entrezgene_manual, $species_id);
	$entrezgene_count++;
	if (defined($array[3])) {     # dead name, add to synonym
	  my @array2 = split(',\s*', $array[3]);
	  foreach my $arr (@array2){
	    XrefParser::BaseParser->add_to_syn($array[0], $hgnc_entrezgene_manual, $arr);
	  }
	}
	
	if (defined($array[4])) {     # alias, add to synonym
	  my @array2 = split(',\s*', $array[4]);
	  foreach my $arr (@array2){
	    XrefParser::BaseParser->add_to_syn($array[0], $hgnc_entrezgene_manual, $arr);
	  }
	}
      }    
    }

    if(defined($array[7])){
      if(defined($entrezgene{$array[7]})){
	$seen = 1;
	XrefParser::BaseParser->add_to_xrefs($entrezgene{$array[7]}, $array[0], '', 
					     $array[1], $array[2], "", $hgnc_entrezgene_mapped, $species_id);
	$entrezgene_count++;
	if (defined($array[3])) {     # dead name, add to synonym
	  my @array2 = split(',\s*', $array[3]);
	  foreach my $arr (@array2){
	    XrefParser::BaseParser->add_to_syn($array[0], $hgnc_entrezgene_mapped, $arr);
	  }
	}
	
	if (defined($array[4])) {     # alias, add to synonym
	  my @array2 = split(',\s*', $array[4]);
	  foreach my $arr (@array2){
	    XrefParser::BaseParser->add_to_syn($array[0], $hgnc_entrezgene_mapped, $arr);
	  }
	}
      }    
    }
    if(!$seen){ # Store to keep descriptions etc
      $self->add_xref($array[0], "", $array[1], $array[2], $source_id, $species_id);      
    }


  } # while HGNC

  $hugo_io->getline();

  print "Loaded a total of " . ($refseq_count + $entrezgene_count) . " HGNC xrefs, $refseq_count from RefSeq curated mappings and $entrezgene_count from EntrezGene mappings\n";

  print "$mismatch xrefs could not be associated via RefSeq or EntrezGene\n";

  return 0; # successful

}

sub rename_url_file{
  return "hugo.txt";
}

1;
    


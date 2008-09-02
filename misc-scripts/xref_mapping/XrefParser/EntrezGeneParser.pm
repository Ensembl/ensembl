package XrefParser::EntrezGeneParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: EntrezGeneParser.pm file <source_id> <species_id>\n\n";
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



  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }
  
  my %species_tax_id = %{$self->get_taxonomy_from_species_id($species_id)};
  

  my $eg_io = $self->get_filehandle($file);
  if ( !defined $eg_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my %seen;


  my $head = $eg_io->getline(); # first record are the headers
  chomp $head;
  my (@arr) = split(/\s+/,$head);
  # process this to the correct indexes to use. (incase they change);

  my $gene_id_index = -2;
  my $gene_symbol_index = -2;
  my $gene_desc_index = -2;
  my $gene_tax_id_index = -2;
  my $gene_synonyms_index = -2;
  foreach (my $i=0; $i<= $#arr; $i++){
    #-1 as first one is "#Format:"q
    if($arr[$i] eq "tax_id"){
      $gene_tax_id_index = $i-1;
    }
    elsif($arr[$i] eq "GeneID"){
      $gene_id_index = $i-1;
    }
    elsif($arr[$i] eq "Symbol"){
      $gene_symbol_index = $i-1;
    }
    elsif($arr[$i] eq "description"){
      $gene_desc_index = $i-1;
    }
    elsif($arr[$i] eq "Synonyms"){
      $gene_synonyms_index = $i-1;
    }
  }
  if( $gene_id_index       == -2 ||
      $gene_symbol_index   == -2 ||
      $gene_desc_index     == -2 ||
      $gene_synonyms_index == -2 ||
      $gene_tax_id_index == -2){
    print "HEADER\n$head\n\n";
    print "Unable to get all the indexes needed\n";
    print "gene_id = $gene_id_index\n";
    print "tax_id = $gene_tax_id_index\n";
    print "symbol = $gene_symbol_index\n";
    print "desc = $gene_desc_index\n";
    print "synonyms = $gene_synonyms_index\n";
    return 0; # this is an error
  }
  my $xref_count = 0;
  my $syn_count  = 0;
  while ( $_ = $eg_io->getline() ) {
    chomp;
    my (@arr) = split(/\t/,$_);
    if(!defined($species_tax_id{$arr[$gene_tax_id_index]})){
      next;
    }
    my $acc    = $arr[$gene_id_index];
    if($seen{$acc}){
      next;
    }
    else{
      $seen{$acc} = 1;
    }
    my $symbol = $arr[$gene_symbol_index];
    my $desc   = $arr[$gene_desc_index];
    $self->add_xref($acc,"",$symbol,$desc,$source_id,$species_id);
    $xref_count++;

    my (@syn) = split(/\|/ ,$arr[$gene_synonyms_index]);
    foreach my $synonym (@syn){
      if($synonym ne "-"){
	$self->add_to_syn($acc, $source_id, $synonym, $species_id);
	$syn_count++;
      }
    }
  }

  $eg_io->close();

  print $xref_count." EntrezGene Xrefs added with $syn_count synonyms\n" if($verbose);
  return 0; #successful
}

 
1;

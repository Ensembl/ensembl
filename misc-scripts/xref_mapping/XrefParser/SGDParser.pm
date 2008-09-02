package XrefParser::SGDParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: SGDParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0]);

}

sub run {

  my $self = shift if (defined(caller(1)));

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
  
 

  my $sgd_io = $self->get_filehandle($file);

  if ( !defined $sgd_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $xref_count =0;
  my $syn_count =0;

  while ( $_ = $sgd_io->getline() ) {
    chomp;
    my ($locus_name, $alias_name, $desc, $gene_prod, $phenotype, $orf_name, $sgd_id) = split(/\t/,$_);

    my (@syn) = split(/\|/,$alias_name);
    $self->add_xref($sgd_id,"",$locus_name,$desc,$source_id,$species_id);
    $xref_count++;
    foreach my $synonym (@syn){
      $self->add_to_syn($sgd_id, $source_id, $synonym, $species_id);
      $syn_count++;
    }
  }

  $sgd_io->close();

  print $xref_count." SGD Xrefs added with $syn_count synonyms\n" if($verbose);
  return 0; #successful
}

1;

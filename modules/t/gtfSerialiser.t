#!/usr/bin/env perl
#
# gtfSerializer.t
#
# Test Bio::EnsEMBL::Utils::IO::GTFSerializer.
#
# Strategy:
#
#   In theory, we should know the GTF output for some 
#   transcripts and compare it against the output of the 
#   serializer. In practice, this is difficult and very 
#   time consuming.
#   One way to circumnvent this is to start from the observation
#   that the serializer was built from the ashes of the web
#   script ensembl-webcode/utils/ensembl2gtf_forweb used by
#   production to get GTFs for last release (71).
#   The strategy is then to download the GTF files for all last release
#   species from the FTP site, and for each file compare its output
#   relative to transcripts with the corresponding output of the 
#   serializer.
#   
# Notes:
# 
#   - This strategy only works for testing the GTFSerializer during
#     the transition from 71 to 72, as after the current release the 
#     GTF files will be produced using the GTFSerializer itself via
#     the Hive pipeline.
#    
#   - The first run of the test reports 9 out 610 failuere for certain
#     species (canis_familiaris cavia_porcellus danio_rerio equus_caballus 
#     loxodonta_africana ornithorhynchus_anatinus oryctolagus_cuniculus 
#     xenopus_tropicalis). This is because the GTFSerializer has corrected
#     a bug of the original script, which was using an incorrect codon
#     table assumption in method check_start_and_stop which is not working
#     for all species.
#     Hence, this failures are ok.
#


use strict;
use warnings;

use Test::More;
use Test::Differences;
use IO::String;
use Net::FTP;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Utils::IO qw/gz_work_with_file/; # iterate_lines
# use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::IO::GTFSerializer;

# my $current_release = software_version();
# my $last_release = $current_release - 1;
my $release = 71;

my $registry = 'Bio::EnsEMBL::Registry';
eval { 
  $registry->load_registry_from_db(
    -host => 'mysql-ensembl-mirror.ebi.ac.uk',
    -port => 4240,
    -db_version => $release,
    -user => 'anonymous',
    -verbose => 0
  )
};
throw "Can't connect to db host" if $@;


SKIP: {
  note "Connecting to Ensembl FTP site";
  my $ftp = eval { ensembl_ftp_connect() };
  skip "Cannot login to Ensembl FTP site. We cannot continue with the tests", 1
    if $@;

  note "Retrieving species list";
  my $species_list = eval { get_species_list($ftp, $release) };
  skip "Cannot retrieve GTF species list from FTP site (release $release). We cannot continue with the tests.", 1
    if $@;
  skip "Empty GTF species list from FTP site (release $release). We cannot continue with the tests", 1
    unless scalar $species_list;

  foreach my $species (@{$species_list}) {
    note "Testing GTF dumps for species $species";
    my $species_gtf = 
      retrieve_gtf_for_species($ftp, $release, $species);  

    skip "Cannot retrieve GTF for $species. We cannot continue with the tests", 1
      unless $species_gtf;

    note "Fetching GTF for 10 randomly chosen transcripts";
    my $transcripts_gtf = eval { get_random_transcripts_gtf($species_gtf, 10) };
    skip "Couldn't fetch 10 transcript. We cannot continue the test.", 1
      unless scalar keys %{$transcripts_gtf} == 10;

    note "Comparing GTF serializer output with release $release data for $species";
    compare_transcript_output($species, $transcripts_gtf);
    
    unlink $species_gtf;
  } 


}

done_testing();

sub compare_transcript_output {
  my ($species, $transcripts_gtf) = @_;

  my $transcript_adaptor =
    $registry->get_adaptor( $species, 'Core', 'Transcript' );

  foreach my $transcript_id (keys %{$transcripts_gtf}) {
    my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);
    defined $transcript or 
      $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);

    skip "Cannot retrieve transcript $transcript_id. Skipping test", 1
      unless defined $transcript;

    my $fh = IO::String->new();
    my $gtf_serializer = 
      Bio::EnsEMBL::Utils::IO::GTFSerializer->new($fh);
    $gtf_serializer->print_feature($transcript);

    is(${$fh->string_ref()}, $transcripts_gtf->{$transcript_id}, "Transcript $transcript_id serialises to GTF as expected");
    # if ($eq_out == 0) {
    #   diag ${$fh->string_ref()};
    #   diag $transcripts_gtf->{$transcript_id};
    # }
    
  }
}

sub get_random_transcripts_gtf {
  my ($gtf_file, $how_many_transcripts) = @_;
  my $transcripts_gtf;

  gz_work_with_file($gtf_file, 'r', sub {
    my ($fh) = @_;

    my $last_insert;
    while( defined $fh and my $line = <$fh> and scalar keys %{$transcripts_gtf} <= $how_many_transcripts) {                              
      $line =~ /transcript_id \"(.+?)\"/;
      $transcripts_gtf->{$1} .= $line;
      $last_insert = $1;
    }
    delete $transcripts_gtf->{$last_insert};
    return;
  });

  return $transcripts_gtf;
}


sub retrieve_gtf_for_species {
  my ($ftp, $release, $species) = @_;

  $ftp->cwd("/pub/release-$release/gtf/$species");
  my @gtf_files = grep { $_ =~ /gtf/ } @{$ftp->ls()};

  throw("No GTF file detected for $species")
    unless scalar @gtf_files;
  throw("More than one GTF file detected for $species")
    if scalar @gtf_files > 1;

  $ftp->binary;
  return $ftp->get($gtf_files[0]);
}

sub get_species_list {
  my ($ftp, $release) = @_;
  $ftp->cwd("/pub/release-$release/gtf") 
    or die "Cannot cd to pub/release-$release/gtf";

  return $ftp->ls();
}

sub ensembl_ftp_connect {
  my $ftp = Net::FTP->new("ftp.ensembl.org", Debug => 0) or
    die "Cannot connect to ensembl.org: $@";
  $ftp->login("anonymous",'-anonymous@') or die "Cannot login ", $ftp->message;

  return $ftp;
}





# my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
# my $dba = $db->get_DBAdaptor('core');

# my $id = 'ENSG00000131044';

# my $ga = $dba->get_GeneAdaptor();
# my $gene = $ga->fetch_by_stable_id($id);

# SKIP: {
#   my $fh = IO::String->new();
#   my $ser = Bio::EnsEMBL::Utils::IO::GTFSerializer->new($fh);
#   # $ser->print_main_header([$gene->feature_Slice()]);
#   $ser->print_feature($gene);
  
#   my $expected = <<'OUT';
# ##gff-version 3
# ##sequence-region   20 30274334 30300924
# OUT
#   #Have to do this outside of the HERETO thanks to tabs
#   $expected .= join("\t", 
#     qw/20  EnsEMBL feature 30274334  30300924  . + ./,
#     'ID=ENSG00000131044;logic_name=ensembl;external_name=C20orf125;description=DJ310O13.1.2 (NOVEL PROTEIN SIMILAR DROSOPHILA PROTEIN CG7474%2C ISOFORM 2 ) (FRAGMENT). [Source:SPTREMBL%3BAcc:Q9BR18];biotype=protein_coding' 
#   );
#   $expected .= "\n";

#   is(${$fh->string_ref()}, $expected, 'Gene serialises to GFF3 as expected');
# }

# done_testing();

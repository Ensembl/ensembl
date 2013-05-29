use strict;
use warnings;

use Test::More;
use Data::Dumper;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Utils::IO qw/gz_work_with_file iterate_lines/;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Utils::IO::GTFSerializer;
use IO::String;
use Net::FTP;

my $current_release = software_version();
my $last_release = $current_release - 1;
my $release = $last_release;

SKIP: {
  my $ftp = eval { ensembl_ftp_connect() };
  skip "Cannot login to Ensembl FTP site. We cannot continue with the tests"
    if $@;

  my $species_list = eval { get_species_list($ftp, $release) };
  skip "Cannot retrieve GTF species list from FTP site (release $release). We cannot continue with the tests." 
    if $@;
  skip "Empty GTF species list from FTP site (release $release). We cannot continue with the tests" 
    unless scalar $species_list;

  foreach my $species (@{$species_list}) {
    my $species_gtf = 
      retrieve_gtf_for_species($ftp, $release, $species);  

    skip "Cannot retrieve GTF for $species. We cannot continue with the tests"
      unless $species_gtf;

    note "Fetching GTF for 10 transcripts randomly chosen";
    my $transcripts_gtf = get_random_transcripts_gtf($species_gtf, 10);
    
    
    unlink $species_gtf;
  } 


}

done_testing();

sub get_random_transcripts_gtf {
  my ($gtf_file, $how_many_transcripts) = @_;
  my $transcripts_gtf;

  gz_work_with_file($gtf_file, 'r', sub {
    my ($fh) = @_;
    while( my $line = <$fh> and scalar keys %{$transcripts_gtf} < $how_many_transcripts) {
      if(random > .5) {
        chomp($line);
        $line =~ /transcript_id \"(.+?)\"/; # double quotes?
  	$transcripts_gtf->{$1} = $line;
      }
    }
  });

  return $transcripts_gtf;
}


sub retrieve_gtf_for_species {
  my ($ftp, $release, $species) = @_;

  $ftp->cwd("/pub/release-$release/gtf/$species");
  my @gtf_files = grep { $_ =~ /gtf/ } @{$ftp->ls()};

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
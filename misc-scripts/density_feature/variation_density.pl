#!/usr/local/ensembl/bin/perl -w

# Calculates the variation density from given core database

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Getopt::Long;

use Data::Dumper;
$Data::Dumper::Maxdepth = 2;

my $max_slices = 100;
my $bin_count  = 150;

my ($host, $user, $pass, $port, $species);

my ($block_count, $genome_size, $block_size );

GetOptions( "host=s",     \$host,
	    "user=s",     \$user,
	    "pass=s",     \$pass,
	    "port=i",     \$port,
	    "species=s",  \$species );

Bio::EnsEMBL::Registry->load_registry_from_db(-host => $host, -user => $user, -pass => $pass, -port => $port);

my $density_feature_adaptor   = Bio::EnsEMBL::Registry->get_adaptor($species, "core", "DensityFeature")        || die "Can't create density feature adaptor";
my $density_type_adaptor      = Bio::EnsEMBL::Registry->get_adaptor($species, "core", "DensityType")           || die "Can't create density type adaptor";
my $analysis_adaptor          = Bio::EnsEMBL::Registry->get_adaptor($species, "core", "analysis")              || die "Can't create analysis adaptor";
my $slice_adaptor             = Bio::EnsEMBL::Registry->get_adaptor($species, "core", "slice")                 || die "Can't create slice adaptor";

my $variation_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, "variation", "VariationFeature") || die "Can't create variation feature adaptor";

# TODO - variation from registry

# Clean up old features first. Also remove analysis and density type entry as these are recreated
my $sth = $slice_adaptor->dbc->prepare("DELETE df, dt, a, ad FROM analysis_description ad, density_feature df, density_type dt, analysis a WHERE ad.analysis_id = a.analysis_id AND a.analysis_id=dt.analysis_id AND dt.density_type_id=df.density_type_id AND a.logic_name='snpdensity'");
$sth->execute();

# Sort slices by coordinate system rank, then by length
my @sorted_slices = sort( {
               $a->coord_system()->rank() <=> $b->coord_system()->rank()
                 || $b->seq_region_length() <=> $a->seq_region_length()
} @{ $slice_adaptor->fetch_all('toplevel') } );

my $analysis =
  new Bio::EnsEMBL::Analysis(
              -program     => "variation_density.pl",
              -database    => "ensembl",
              -gff_source  => "variation_density.pl",
              -gff_feature => "density",
              -logic_name  => "snpdensity",
              -description => 'Density of Single Nucleotide Polymorphisms (SNPs) calculated by variation_density.pl (see scripts at the <a rel="external" href="http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/?root=ensembl">Sanger Institute CVS</a> repository).',
              -display_label => 'SNP Density',
              -displayable   => 1 );

$analysis_adaptor->store($analysis);
$analysis_adaptor->update($analysis);

# Create and store new density type

my $dt = Bio::EnsEMBL::DensityType->new(-analysis        => $analysis,
					-region_features => $bin_count,
					-value_type      => 'sum');
$density_type_adaptor->store($dt);

# Now the actual feature calculation loop

my $slice_count = 0;

my ($current, $current_start, $current_end);

# prepare statement outside of loop
$sth = $variation_feature_adaptor->prepare("SELECT COUNT(*) FROM variation_feature WHERE seq_region_id = ? AND seq_region_start < ? AND seq_region_end > ?");

foreach my $slice (@sorted_slices){

  $block_size = $slice->length() / $bin_count;

  print "Calculating SNP densities for ". $slice->seq_region_name() . " with block size $block_size\n";

  $current_end = 0;
  $current = 0;

  while($current_end < $slice->length) {

    $current += $block_size;
    $current_start = $current_end+1;
    $current_end = int( $current + 1 );

    if ($current_end < $current_start) {
      $current_end = $current_start;
    }

    if ($current_end > $slice->end()) {
      $current_end = $slice->end();
    }

    # Get count of SNPs in this region
    $sth->execute($slice->get_seq_region_id(), $current_end, $current_start);
    my $count = ($sth->fetchrow_array())[0];

    my $df = Bio::EnsEMBL::DensityFeature->new(-seq_region    => $slice,
					       -start         => $current_start,
        				       -end           => $current_end,
        				       -density_type  => $dt,
        				       -density_value => $count);
    $density_feature_adaptor->store($df);

  }

  last if ( $slice_count++ > $max_slices );

}


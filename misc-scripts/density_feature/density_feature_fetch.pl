#!/usr/local/ensembl/bin/perl -w

#
# Ignore file anme this should calculate the gene count for a chromosome (X)
#


use strict;
use lib '/tmp_mnt/nfs6/vol_vol1_homes/ianl/CVS/ensembl/modules';
use lib '/tmp_mnt/nfs6/vol_vol1_homes/ianl/read_only/bioperl-live';


use Data::Dumper;

#use dbi;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DensityTypeAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
print "hello\n";

my $host = '127.0.0.1';
my $user = 'ensro';
#my $pass = 'ensembl';
my $dbname = 'homo_sapiens_core_20_34';
my $port = '5050';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
#					    -pass => $pass,
					    -dbname => $dbname);


my $slice_adaptor = $db->get_SliceAdaptor();
my $dfa = $db->get_DensityFeatureAdaptor();


my $num_of_blocks = 20;
my $interpolate = 1;


#comment out as appropriate.
#my $type="GeneDensity";
#my $type="PercentGC";
my $type ="PercentageRepeat";

foreach my $chrom (qw(X Y)){
  if($interpolate){
    print "getting density feature for chromosome $chrom with $num_of_blocks blocks\n";
  }
  else{
    print "getting ALL  density feature for chromosome $chrom\n";
  }
  my $slice = $slice_adaptor->fetch_by_region('chromosome',$chrom);


  my @features = @{$dfa->fetch_all_by_Slice($slice,$type, $num_of_blocks, $interpolate)};

  print_features(\@features);
}


#
# helper to draw an ascii representation of the density features
#
sub print_features {
  my $features = shift;

  return if(!@$features);

  my $sum = 0;
  my $length = 0;
#  my $type = $features->[0]->{'density_type'}->value_type();

  print("\n");
  my $max=0;
  foreach my $f (@$features) {
    if($f->density_value() > $max){
      $max=$f->density_value();
    }
  }
  foreach my $f (@$features) {
    my $i=1;
    for(; $i< ($f->density_value()/$max)*40; $i++){
      print "*";
    }
    for(my $j=$i;$j<40;$j++){
      print " ";
    }
    print "  ".int($f->density_value())."\t".$f->start()."\t".$f->end()."\n";
  }
}

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;

my $host = 'ecs3d';
my $user = 'ensadmin';
my $pass = 'ensembl';
my $dbname = 'mcvicker_homo_sapiens_core_20_34b';
my $port = '3307';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);


#
# choose a block size that makes 150 blocks for the shortest chromosome
#
my $block_size;
my $min_chr_len;
my $slice_adaptor = $db->get_SliceAdaptor();
my @chromosomes = @{$slice_adaptor->fetch_all('chromosome')};

#my @chromosomes = ($slice_adaptor->fetch_by_region('chromosome', '1'));

foreach my $chr (@chromosomes) {
  next if($chr->seq_region_name =~ /DR/i);
  $min_chr_len = $chr->seq_region_length() if(!defined($min_chr_len) || $chr->seq_region_length < $min_chr_len);
}

$block_size = int($min_chr_len / 150);


#
# Get the adaptors needed;
#

my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();

foreach my $known (1, 0) {
  #
  # Create new analysis object for density calculation.
  #

  my $analysis;

  if($known) {
    $analysis = $aa->fetch_by_logic_name('kngene');
  } else {
    $analysis = $aa->fetch_by_logic_name('gene');
  }

  #
  # Create new density type.
  #

  my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					  -block_size => $block_size,
					  -value_type => 'sum');

  $dta->store($dt);


  foreach my $slice (@chromosomes){
    print "creating density feature for chromosome ".$slice->seq_region_name()."with block size of $block_size\n";

    my $start = $slice->start();
    my $end = ($start + $block_size)-1;
    my $term = $slice->start+$slice->length;

    my @density_features=();
    while($start < $term){

      my $sub_slice = $slice_adaptor->fetch_by_region('chromosome',$slice->seq_region_name(),$start,$end);

      my $count =0;

      #
      # Store info for genes (ignore pseudo genes)
      #

      foreach my $gene (@{$sub_slice->get_all_Genes()}){
	if($gene->analysis()->logic_name() ne "pseudogene" and $gene->start >=1 ){
	  $count++ if(!$known || $gene->is_known());
	}
      }

      push @density_features, Bio::EnsEMBL::DensityFeature->new(-seq_region    => $slice,
								-start         => $start,
								-end           => $end,
								-density_type  => $dt,
								-density_value => $count);

      $start = $end+1;
      $end   = ($start + $block_size)-1;
    }
    $dfa->store(@density_features);
    print scalar @density_features;
    print_features(\@density_features);
  }
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
    print "  ".$f->density_value()."\t".$f->start()."\n";
  }
#  my $avg = undef;
#  $avg = $sum/$length if($length < 0);
#  print("Sum=$sum, Length=$length, Avg/Base=$sum");
}




  



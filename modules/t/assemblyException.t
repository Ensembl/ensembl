use strict;

BEGIN { $| = 1;
	use Test;
	plan tests => 7;
}


use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok(1);

my $db = $multi->get_DBAdaptor('core');
my $sfa = $db->get_SimpleFeatureAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

##chromosome Y is a fake 'PAR' linked to chromosome 20
my $slice = $slice_adaptor->fetch_by_region('chromosome', 'Y',8e6,13e6);
my $feats = $sfa->fetch_all_by_Slice($slice);
debug("Got " . scalar(@$feats));
ok( @$feats ==58 );

print_features($feats);


#HAP_1 is a fake haplotype on chromosome 20
$slice = $slice_adaptor->fetch_by_region('chromosome', '20_HAP1',
                                            30_399_998,30_600_002);
my $org_slice = $slice_adaptor->fetch_by_region('chromosome', '20',
                                            30_430_000,30_500_000 );


$feats = $sfa->fetch_all_by_Slice($slice);

debug("Got " . scalar(@$feats));
ok( @$feats == 9 );

print_features($feats);

$multi->hide( "core", "simple_feature" );
$multi->save( "core", "meta_coord" );

for my $f ( @$feats ) {
  $f->dbID( undef );
  $f->adaptor( undef );
  $sfa->store( $f );

  $f->dbID( undef );
  $f->adaptor( undef );
  $f->slice( $org_slice );
  $sfa->store( $f );
    
}


$slice = $slice_adaptor->fetch_by_region('chromosome', '20_HAP1',
					 30_400_000,30_600_000);
$feats = $sfa->fetch_all_by_Slice( $slice );

debug( "After storing retrieval" );
print_features($feats);
ok(@$feats == 14);



#
# sequence getting tests
#

my $hap_slice = $slice_adaptor->fetch_by_region('chromosome', '20_HAP1',
                                             30_400_000,30_700_000 );

$org_slice = $slice_adaptor->fetch_by_region('chromosome', '20',
                                             30_400_000,30_800_000 );

my ( $fhs, $bhs, $fos, $bos );

debug( "Front hap seq: ".($fhs = $hap_slice->subseq( 99_991, 100_000 )));
debug( "Back hap seq: ".($bhs = $hap_slice->subseq( 200_001, 200_010 )));
debug( "Front org seq: ".( $fos = $org_slice->subseq( 99_991, 100_000 )));
debug( "Back org seq: ".( $bos = $org_slice->subseq( 300_001, 300_010 )));

ok( $fhs eq $fos );
ok( $bhs eq $bos );


$slice = $slice_adaptor->fetch_by_region('chromosome', '20_HAP1',
					 30_499_998,30_500_002);

debug($slice->seq);
ok( $slice->seq() eq "GTNNN" );

$multi->restore();


#try projecting a hap slice to the contig coordinate system
debug("Org slice projection");
my $projection = $org_slice->project('contig');
print_projection($projection);

debug("Hap slice projection");
$projection = $hap_slice->project('contig');
print_projection($projection);


sub print_features {
  return if(!$verbose);
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      my $seqname = $f->slice->seq_region_name();
      my $analysis = $f->analysis->logic_name();
      debug($seqname . ' ' . $f->start().'-'.$f->end().'('.$f->strand().
            ') ['. $f->dbID.'] '.$f->score() ." ($analysis) ");
    } else {
      debug('UNDEF');
    }
  }
}



sub print_projection {
  my $proj = shift;
  foreach my $seg (@$proj) {
    my ($start, $end, $seq_reg) = ($seg->[0],$seg->[1],$seg->[2]->seq_region_name());
    debug("[$start-$end] $seq_reg");
  }
}

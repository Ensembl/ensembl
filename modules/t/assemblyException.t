use lib 't';
use strict;

BEGIN { $| = 1;
	use Test;
	plan tests => 5;
}


use MultiTestDB;

use Bio::EnsEMBL::DnaDnaAlignFeature;
use TestUtils qw(test_getter_setter debug);

our $verbose = 0;

my $multi = MultiTestDB->new();
ok(1);

my $db = $multi->get_DBAdaptor('core');
my $sfa = $db->get_SimpleFeatureAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

##chromosome Y is a fake 'PAR' linked to chromosome 20
my $slice = $slice_adaptor->fetch_by_region('chromosome', 'Y',8e6,13e6);
my $feats = $sfa->fetch_all_by_Slice($slice);
debug("Got " . scalar(@$feats));
print_features($feats);


#HAP_1 is a fake haplotype on chromosome 20
my $slice = $slice_adaptor->fetch_by_region('chromosome', '20_HAP1',
                                            30_399_998,30_600_002);
my $org_slice = $slice_adaptor->fetch_by_region('chromosome', '20',
                                            30_430_000,30_500_000 );


my $feats = $sfa->fetch_all_by_Slice($slice);

debug("Got " . scalar(@$feats));
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

my $org_slice = $slice_adaptor->fetch_by_region('chromosome', '20',
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

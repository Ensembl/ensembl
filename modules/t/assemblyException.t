use lib 't';
use strict;

BEGIN { $| = 1;
	use Test;
	plan tests => 4;
}


use MultiTestDB;

use Bio::EnsEMBL::DnaDnaAlignFeature;
use TestUtils qw(test_getter_setter debug);

our $verbose = 0;

my $multi = MultiTestDB->new();
ok(1);

my $db = $multi->get_DBAdaptor('core');
my $dafa = $db->get_DnaAlignFeatureAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

##chromosome Y is a fake 'PAR' linked to chromosome 20
#my $slice = $slice_adaptor->fetch_by_region('chromosome', 'Y',8e6,13e6);
#my $feats = $dafa->fetch_all_by_Slice($slice);
#debug("Got " . scalar(@$feats));
#print_features($feats);


#HAP_1 is a fake haplotype on chromosome 20
my $slice = $slice_adaptor->fetch_by_region('chromosome', '20_HAP1',
                                            30_499_998,30_500_002);

my $feats = $dafa->fetch_all_by_Slice($slice);
debug("Got " . scalar(@$feats));
print_features($feats);

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


debug($slice->seq);
ok( $slice->seq() eq "GTNNN" );


sub print_features {
  return if(!$verbose);
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      my $seqname = $f->slice->seq_region_name();
      my $analysis = $f->analysis->logic_name();
      debug($seqname . ' ' . $f->start().'-'.$f->end().'('.$f->strand().
            ') ['. $f->dbID.'] '. $f->cigar_string . ' ' .
            $f->hstart .'-'.$f->hend.' ('.$f->hstrand.')'.$f->score() .
            " ($analysis) " . $f->percent_id);
    } else {
      debug('UNDEF');
    }
  }
}

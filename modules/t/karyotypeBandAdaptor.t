use strict;
use warnings;


BEGIN { $| = 1;
	use Test;
	plan tests => 18;
}


use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::KaryotypeBand;

our $verbose= 0;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');


my $kba = $db->get_KaryotypeBandAdaptor();

ok(ref($kba) && $kba->isa('Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor'));


my $slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                   '20',
                                                   30_000_000,
                                                   32_000_000);




my $bands = $kba->fetch_all_by_Slice($slice);

print_features($bands);

ok(@$bands == 1);

my ($band) = @$bands;

ok($band->name() eq 'q11.21');
ok($band->start() == 200001);
ok($band->end()   == 2600001);
ok($band->strand() == 0);



($band) = @{$kba->fetch_all_by_chr_name('20')};

print_features([$band]);

ok($band->name() eq 'q11.21');
ok($band->slice->seq_region_name eq '20');
ok($band->start() == 30_200_000);
ok($band->end()   == 32_600_000);


($band) = @{$kba->fetch_all_by_chr_band('20', 'q11')};

print_features([$band]);

ok($band->name() eq 'q11.21');
ok($band->slice->seq_region_name eq '20');
ok($band->start() == 30_200_000);
ok($band->end()   == 32_600_000);

$band = $kba->fetch_by_dbID(1);
print_features([$band]);
ok($band->name() eq 'q11.21');
ok($band->slice->seq_region_name eq '20');
ok($band->start() == 30_200_000);
ok($band->end()   == 32_600_000);

sub print_features {
  my $features = shift;

  foreach my $f (@$features) {
    if(defined($f)) {
      debug('['.$f->dbID.']'.$f->start.'-'.$f->end.'('.$f->strand.')'.
            '  '.$f->name(). ' '.$f->stain());
    } else {
      debug('UNDEF');
    }
  }
}

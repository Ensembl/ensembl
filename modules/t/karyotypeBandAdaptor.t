# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Test::Warnings;

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


# Need complete slice to store band as coordinates will be relative to that slice
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', '20');
my $new_band = Bio::EnsEMBL::KaryotypeBand->new
  (-START       => 30_000_000,
   -END         => 30_200_000,
   -SLICE       => $slice,
   -NAME        => 'test_band');

$kba->store(($new_band));

is(scalar(@{$kba->fetch_all_by_chr_name('20')}), 2, 'Added one more band');

$new_band = @{$kba->fetch_all_by_chr_name('20')}[1];

is($new_band->name(),  'test_band', 'New band name matches');
is($new_band->start(), 30_000_000, 'New band start matches');
is($new_band->end(), 30_200_000, 'New band end matches');
is($new_band->strand(), 0, 'Bands have no strand');



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

done_testing();

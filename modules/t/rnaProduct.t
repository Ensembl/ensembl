# Copyright [2018] EMBL-European Bioinformatics Institute
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

## no critic (RequireFilenameMatchesPackage)

package RNAProductTests;

use strict;
use warnings;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::RNAProduct;
use Bio::EnsEMBL::Transcript;

use Test::More;
use Test::Warnings;
use Test::Exception;

my $loaded = 0;
END { print "not ok 1 - Test set-up completed\n" unless $loaded; }

use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

$loaded = 1;

ok(1, 'Test set-up completed');

my $db = $multi->get_DBAdaptor('core');

my $rp = Bio::EnsEMBL::RNAProduct->new();

ok($rp, 'RNAProduct constructor works without arguments');


# We will use the minimally constructed object from above for further
# testing so let us get rid of this one as soon as we are done with it.
{
  my %cta = (
    start => 123,
    end => 456,
    stable_id => 'ENSM00012345',
    version => 1337,
    dbID => 314,
    seq => 'ACGTACGT',
    created_date => time(),
    modified_date => time()
  );
  my $rp_with_args = Bio::EnsEMBL::RNAProduct->new(
    -SEQ_START => $cta{start},
    -SEQ_END => $cta{end},
    -STABLE_ID => $cta{stable_id},
    -VERSION => $cta{version},
    -DBID => $cta{dbID},
    -SEQ => $cta{seq},
    -CREATED_DATE => $cta{created_date},
    -MODIFIED_DATE => $cta{modified_date}
  );
  foreach my $member (sort keys %cta) {
    is($rp_with_args->{$member}, $cta{$member}, "RNAProduct constructor sets $member correctly");
  }
}

is($rp->version(), 1, 'Default rnaproduct version == 1');

ok(test_getter_setter($rp, 'start', 42), 'Test getter/setter start()');
ok(test_getter_setter($rp, 'end', 64), 'Test getter/setter end()');
ok(test_getter_setter($rp, 'stable_id', 1), 'Test getter/setter stable_id()');
ok(test_getter_setter($rp, 'dbID', 3), 'Test getter/setter dbID()');
ok(test_getter_setter($rp, 'version', 13), 'Test getter/setter version()');
ok(test_getter_setter($rp, 'created_date', time()), 'Test getter/setter created_date()');
ok(test_getter_setter($rp, 'modified_date', time()), 'Test getter/setter modified_date()');
ok(test_getter_setter($rp, 'seq', 'AATTGGCC'), 'Test getter/setter seq()');

subtest 'Test stable_id_version() functionality' =>  sub {
  ok(test_getter_setter($rp, 'stable_id_version', 3.14),
     'getter/setter with \'stable_id.version\' as input');
  ok(test_getter_setter($rp, 'stable_id_version', 'aqq'),
     'getter/setter with \'stable_id\' as input');

  # Let's be paranoid and assume test_getter_setter() cleans up after itself
  $rp->stable_id_version('ENSfoo.1');
  is($rp->stable_id_version(), $rp->stable_id() . '.' . $rp->version(),
     'set by stable_id_version(), get by stable_id() + version()');

  $rp->stable_id('ENSbar');
  $rp->version(9);
  is($rp->stable_id_version(), $rp->stable_id() . '.' . $rp->version(),
     'set by stable_id() + version(), get by stable_id_version()');
};

subtest 'display_id() functionality' =>  sub {
  # Start with a minimal object and gradually add missing data
  my $rp_blank = Bio::EnsEMBL::RNAProduct->new();

  is($rp_blank->display_id(), '',
     'return empty string if neither stable_id nor dbID exist');

  $rp_blank->dbID(12345);
  is($rp_blank->display_id(), $rp_blank->dbID(),
     'return dbID if no stable_id exists');

  $rp_blank->stable_id(54321);
  is($rp_blank->display_id(), $rp_blank->stable_id(),
     'return stable_id if it exists');
};

{
  # Again, assume test_getter_setter() cleans up after itself
  my $dummy_sequence = 'CGATCCGGAAAA';
  $rp->seq($dummy_sequence);
  is($rp->length(), length($dummy_sequence), 'Check if length() returns correct value');
}

{
  dies_ok(sub { $rp->transcript({ }) }, 'Transcript setter dies on incorrect argument type');
  # The first step of test_getter_setter() is to preserve the existing value
  # of the property being tested. This is normally fine but in case of
  # transcript() invoking it as getter with no transcript previously having
  # been set causes it to attempt a database lookup, which:
  #  - is not desired because we are testing local functionality now, and
  #  - will fail because no adaptor has been set either.
  # Therefore, set a dummy transcript first. Then use a *different* dummy
  # transcript in the test to make sure the setter really works.
  # Nb. it is necessary to assign the first dummy transcript to a variable
  # because the setter uses a weak reference. If you simply use the return
  # value of the Transcript constructor as an argument to transcript(), the
  # former will get garbage-collected and the transcript will remain unset.
  my $dummy_transcript = Bio::EnsEMBL::Transcript->new();
  $rp->transcript($dummy_transcript);
  ok(test_getter_setter($rp, 'transcript', Bio::EnsEMBL::Transcript->new()), 'Test getter/setter transcript()');
}

is($rp->cdna_start(), $rp->start(),
   'Test if cdna_start() returns the same value as start()');
is($rp->cdna_end(), $rp->end(),
   'Test if cdna_end() returns the same value as end()');


# TODO: More RNAProduct tests


#
# Tests for the mature-RNA adaptor
##################################

my $rp_a  = $db->get_RNAProductAdaptor();

ok($rp_a, 'Can get RNAProductAdaptor from core DBAdaptor');

is($rp_a->fetch_by_dbID(0), undef, 'Adaptor returns undef for nonexistent dbID');
is($rp_a->fetch_by_stable_id('fnord'), undef, 'Adaptor returns undef for nonexistent stable ID');

subtest 'fetch_all_by_Transcript() functionality' => sub {
  my $t_a = $db->get_TranscriptAdaptor();
  my $rps;
  my $t;

  $t = Bio::EnsEMBL::Transcript->new();
  $rps = $rp_a->fetch_all_by_Transcript($t);
  isa_ok($rps, 'ARRAY', 'fetch_all_by_Transcript() return value');
  is(scalar @{$rps}, 0, 'Adaptor returns empty list for invalid Transcript');

  $rps = undef;
  $t = $t_a->fetch_by_dbID(21716);
  $rps = $rp_a->fetch_all_by_Transcript($t);
  is(scalar @{$rps}, 0, 'Adaptor returns empty list for Transcript with no RNA products');

  $rps = undef;
  $t = $t_a->fetch_by_dbID(21717);
  $rps = $rp_a->fetch_all_by_Transcript($t);
  # Do not bother checking if elements of the returned array are defined,
  # we are testing the method and not database consistency.
  cmp_ok(scalar @{$rps}, '>', 0, 'Non-empty list for Transcript with RNA products');
};

$rp = undef;
$rp = $rp_a->fetch_by_dbID(1);
ok($rp, 'Can fetch RNAProduct by dbID');

$rp = undef;
$rp = $rp_a->fetch_by_stable_id('ENSM00000000001');
ok($rp, 'Can fetch RNAProduct by stable ID');

# FIXME: perform an in-depth inspection of one of the fetched RNAProducts,
# to make sure new_fast() call all of these fetch methods use does what it
# is supposed to do.


# TODO: More RNAProductAdaptor tests


# Test generic_count(), inherited method from BaseAdaptor
is($rp_a->generic_count(), @{$rp_a->list_dbIDs()}, "Number of features from generic_count is equal to the number of dbIDs from list_dbIDs");

done_testing();

1;

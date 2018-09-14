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
use Bio::EnsEMBL::MicroRNA;
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
my $type_mapper = Bio::EnsEMBL::Utils::RNAProductTypeMapper::mapper();


#
# Tests for the RNAProductTypeMapper
#

subtest 'RNAProductTypeMapper tests' => sub {
  my $rpt_mapper = Bio::EnsEMBL::Utils::RNAProductTypeMapper->mapper();

  my $rpt_mapper2 = Bio::EnsEMBL::Utils::RNAProductTypeMapper->mapper();
  is($rpt_mapper, $rpt_mapper2, 'mapper() reuses existing instance if present');

  is($rpt_mapper->type_code_to_class('miRNA'), 'Bio::EnsEMBL::MicroRNA',
     'Can map existing type ID to class');
  dies_ok(sub { $rpt_mapper->type_code_to_class('semprini'); },
	  'Exception thrown on unknown type ID');

  is($rpt_mapper->class_to_type_code('Bio::EnsEMBL::RNAProduct'), 'generic',
     'Can map existing class to type ID');
  dies_ok(sub { $rpt_mapper->class_to_type_code('Bio::EnsEMBL::Storable'); },
	  'Exception thrown on unknown rnaproduct class name');
};


#
# Tests for offline RNAProduct objects
#

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

is($rp->type_code(), $type_mapper->class_to_type_code(ref($rp)),
   'RNAProduct object has expected type code');

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
  dies_ok(sub { $rp->seq() }, 'Sequence getter dies if neither local nor DB data is available');
  # Again, assume test_getter_setter() cleans up after itself
  my $dummy_sequence = 'CGATCCGGAAAA';
  $rp->seq($dummy_sequence);
  is($rp->length(), length($dummy_sequence), 'Check if length() returns correct value');
  ok(test_getter_setter($rp, 'seq', 'AACCGGTT'), 'Test getter/setter seq()');
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


#
# Tests for the mature-RNA adaptor
#

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

subtest 'fetch_all_by_type() functionality' => sub {
  my $n_rps;

  # At the moment we have only got miRNA in the homo_sapiens test database

  # FIXME: compare this to the total number of RNAProducts?
  $n_rps = scalar @{$rp_a->fetch_all_by_type('miRNA')};
  cmp_ok($n_rps, '>', 0, 'Got non-empty list of miRNA rnaproducts');
  $n_rps = scalar @{$rp_a->fetch_all_by_type('generic')};
  cmp_ok($n_rps, '==', 0, 'Got empty list of generic rnaproducts');
};

$rp = undef;
$rp = $rp_a->fetch_by_dbID(1);
ok($rp, 'Can fetch RNAProduct by dbID');

$rp = undef;
$rp = $rp_a->fetch_by_stable_id('ENSM00000000001');
ok($rp, 'Can fetch RNAProduct by stable ID');

isa_ok($rp, 'Bio::EnsEMBL::MicroRNA', 'miRNA object from database');
is($rp->type_code(), $type_mapper->class_to_type_code(ref($rp)),
   'MicroRNA object has expected type code');

# FIXME: perform an in-depth inspection of one of the fetched RNAProducts,
# to make sure new_fast() call all of these fetch methods use does what it
# is supposed to do.

is($rp->seq(), 'AAAAACCCAGGAATCACCTGGA', 'Can retrieve associated sequence');

# Do not check any data inside the Transcript object, it is not our job to
# check database consistency. Just check that we do get something back.
isnt($rp->transcript(), undef, 'Can retrieve associated Transcript object');

# And now, force the Transcript association to be built on the fly.
$rp->transcript(undef);
isnt($rp->transcript(), undef, 'Transcript association can be built on demand for valid dbID');

# FIXME: might want to add tests for the reverse strand as well
is($rp->genomic_start(), $rp->transcript()->start() + $rp->start() - 1,
   'genomic_start() gives correct values (forward strand)');
is($rp->genomic_end(), $rp->transcript()->start() + $rp->end() - 1,
   'genomic_end() gives correct values (forward strand)');


subtest 'Attribute functionality' => sub {
  my $rp_all_attrs = $rp->get_all_Attributes();
  cmp_ok(scalar @$rp_all_attrs, '>', 0, 'Get a non-empty list of attributes');

  my $rp_notes = $rp->get_all_Attributes('note');
  cmp_ok(scalar @$rp_notes, '>', 0, 'Get a non-empty list of \'note\' attributes');

  my $rp_nonsense = $rp->get_all_Attributes('xyzzy');
  is(scalar @$rp_nonsense, 0, 'Get empty attribute list for nonsense code');

  dies_ok(sub { $rp->add_Attributes({}) },
	  'add_Attributes() dies on invalid argument type');

  my $n_attrs_before = scalar @$rp_all_attrs;
  my $extra_attr1 = Bio::EnsEMBL::Attribute->new(
    -CODE => 'note',
    -NAME => 'Note',
    -VALUE => 'and another thing'
  );
  my $extra_attr2 = Bio::EnsEMBL::Attribute->new(
    -CODE => '_rna_edit',
    -VALUE => '1 6 GATTACA',
    -NAME => 'RNA editing'
  );
  $rp->add_Attributes($extra_attr1, $extra_attr2);
  is(scalar @{$rp->get_all_Attributes()}, scalar $n_attrs_before + 2, 'Added two new attributes');

  # FIXME: Add SeqEdit tests once we have got some meaningful data for this
  # in the test database. The way this is done in Transcript tests ought to
  # be a good reference.
};


subtest 'xref functionality' => sub {
  my $xrefs = $rp->get_all_DBEntries();
  cmp_ok(scalar @$xrefs, '>', 0, 'Got a non-empty list of DBEntries');

  dies_ok(sub { $rp->add_DBEntry({}) },
	  'add_DBEntry() dies on invalid argument type');

  my $n_xrefs_before = scalar @$xrefs;
  my $dbe = Bio::EnsEMBL::DBEntry->new(
    -primary_id => 'test_id',
    -version    => 1,
    -dbname     => 'miRNA_Registry',
    -display_id => 'test_id'
  );
  $rp->add_DBEntry($dbe);
  is(scalar @{$rp->get_all_DBEntries()}, scalar $n_xrefs_before + 1, 'Added one new xref');

  # No need for deep comparisons here, these four methods are supposed
  # to return literally the same reference
  is($xrefs, $rp->get_all_object_xrefs(),
     'get_all_object_xrefs() is a alias of get_all_DBEntries()');
  is($xrefs, $rp->get_all_DBLinks(),
     'get_all_DBLinks() is a alias of get_all_DBEntries()');
  is($xrefs, $rp->get_all_xrefs(),
     'get_all_xrefs() is a alias of get_all_DBEntries()');
};

my $rp_exts = $rp_a->fetch_all_by_external_name('hsa-miR-1-3p');
cmp_ok(scalar @$rp_exts, '>', 0, 'Can fetch RNAProduct by external ID');

# Test generic_count(), inherited method from BaseAdaptor
is($rp_a->generic_count(), @{$rp_a->list_dbIDs()}, "Number of features from generic_count is equal to the number of dbIDs from list_dbIDs");


#
# MicroRNA-specific tests
#

subtest 'MicroRNA tests' => sub {
  my $mirna;

  $mirna = Bio::EnsEMBL::MicroRNA->new();
  ok($mirna, 'MicroRNA constructor works without arguments');
  isa_ok($mirna, 'Bio::EnsEMBL::MicroRNA', 'miRNA object from new()');

  $mirna = Bio::EnsEMBL::MicroRNA->new(
    -SEQ_START => 314,
    -ARM => 3
  );
  ok($mirna, 'MicroRNA constructor works with arguments');
  is($mirna->arm(), 3, 'MicroRNA-specific parameters set OK');
  is($mirna->start(), 314, 'Generic RNAProduct parameters set OK');

  $mirna = $rp_a->fetch_by_dbID(1);
  ok($mirna, 'Can fetch MicroRNA from RNAProductAdaptor');
  isa_ok($mirna, 'Bio::EnsEMBL::MicroRNA', 'miRNA object from RNAProductAdaptor');
  isnt($mirna->arm(), undef, 'Can retrieve miRNA arm value from DB');
};


#
# Operations which modify database entries
#

subtest 'Write operations' => sub {
  my $rp_count = count_rows($db, 'rnaproduct');
  my $rp_attr_count = count_rows($db, 'rnaproduct_attrib');
  my $object_xref_count = count_rows($db, 'object_xref');
  my $xref_count = count_rows($db, 'xref');
  $multi->save('core', 'rnaproduct', 'rnaproduct_attrib', 'object_xref', 'xref');

  my $t_a = $db->get_TranscriptAdaptor();
  my $ins_parent = $t_a->fetch_by_dbID(21728);

  my %insert_args = (
    start         => 6,
    end           => 27,
    stable_id     => 'ENSM00000000099',
    version       => 1,
    created_date  => time(),
    modified_date => time(),
    arm           => 5,
  );
  my $ins_rp = Bio::EnsEMBL::MicroRNA->new(
    -SEQ_START     => $insert_args{start},
    -SEQ_END       => $insert_args{end},
    -STABLE_ID     => $insert_args{stable_id},
    -VERSION       => $insert_args{version},
    -CREATED_DATE  => $insert_args{created_date},
    -MODIFIED_DATE => $insert_args{modified_date},
    -ARM           => $insert_args{arm},
  );
  my $ins_dbe = Bio::EnsEMBL::DBEntry->new(
    -primary_id => 'MIMAT0000062',
    -version    => 1,
    -dbname     => 'miRNA_Registry',
    -display_id => 'hsa-let-7a-5p',
  );
  is($ins_rp->dbID(), undef, 'MicroRNA dbID not defined before storing');
  $ins_rp->add_DBEntry($ins_dbe);
  $rp_a->store($ins_rp, $ins_parent->dbID());
  isnt($ins_rp->dbID(), undef, 'MicroRNA dbID defined after storing');
  is($ins_rp->adaptor(), $rp_a,
     'MicroRNA linked to RNAPRoductAdaptor after storing');
  is(count_rows($db, 'rnaproduct'), $rp_count + 1,
     'Post-store rnaproduct row count correct');
  is(count_rows($db, 'rnaproduct_attrib'), $rp_attr_count + 1,
     'Post-store rnaproduct_attrib row count correct');
  is(count_rows($db, 'object_xref'), $object_xref_count + 1,
     'Post-store object_xref row count correct');
  is(count_rows($db, 'xref'), $xref_count + 1,
     'Post-store xref row count correct');

  my $fetched_rp = $rp_a->fetch_by_dbID($ins_rp->dbID());
  for my $method (keys %insert_args) {
    is($fetched_rp->$method(), $insert_args{$method},
       "RNAProduct from database has correct $method value");
  }
  is(scalar @{ $fetched_rp->get_all_Attributes() }, 1,
     'Fetched MicroRNA has correct number of attributes');
  is(scalar @{ $fetched_rp->get_all_DBEntries() }, 1,
     'Fetched MicroRNA has correct number of xrefs');
  isnt($fetched_rp->transcript(), undef,
       'Fetched MicroRNA has a parent transcript');

  my $del_rp = $rp_a->fetch_by_dbID(1);
  $rp_a->remove($del_rp);
  is($del_rp->dbID(), undef, 'Removed MicroRNA has no dbID');
  is($del_rp->adaptor(), undef, 'Removed MicroRNA not linked to adaptor');
  is(count_rows($db, 'rnaproduct'), $rp_count,
     'Post-remove rnaproduct row count correct');
  is(count_rows($db, 'rnaproduct_attrib'), $rp_attr_count - 1,
     'Post-remove rnaproduct_attrib row count correct');
  is(count_rows($db, 'object_xref'), $object_xref_count,
     'Post-remove object_xref row count correct');
  # Not a mistake, actual xrefs do NOT get removed by
  # DBEntryAdaptor::remove_from_object(). See e.g. Doxygen docs
  # for that method.
  is(count_rows($db, 'xref'), $xref_count + 1,
     'Post-remove xref row count correct');
  is ($rp_a->fetch_by_dbID(1), undef,
      'Fetching removed MicroRNA returns undef');

  $multi->restore('core', 'rnaproduct', 'rnaproduct_attrib', 'object_xref', 'xref');
};


done_testing();

1;

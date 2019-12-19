=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Xref::Test::TestDB;
use FindBin '$Bin';

my $db = Xref::Test::TestDB->new();
ok($db, 'default instantiation proceeds as planned with testdb.conf');
throws_ok { Xref::Test::TestDB->new(config_file => 'not_here') }
  qr/does not exist!/,
  'TestDB grumbles differently when an invalid config file is offered';

# This auto-deploys the schema
$db = Xref::Test::TestDB->new(
  # config_file => 'testdb.conf'
  config => {
    driver => 'SQLite',
    file => 'test.db',
    create => 1
  }
);

ok($db, 'TestDB ready to go');

my $source = $db->schema->resultset('Source')->create({
  name => 'RefSeq',
  source_release => '38',
  priority => 1,
  priority_description => 'Like a boss',
});

ok(defined $source->source_id, 'Was the source created in the DB?');


my $xref = $source->create_related('xrefs', {
  accession => 'NM01234',
  version => 1,
  label => 'NM01234.1',
  description => 'Fake RefSeq transcript',
  species_id => '9606',
  info_type => 'DIRECT',
  info_text => 'These are normally aligned',
  dumped => 'NO_DUMP_ANOTHER_PRIORITY'
});

my $rs = $db->schema->resultset('Xref')->search(
  { accession => 'NM01234'}
);

my $matching_xref = $rs->next;
ok(defined $matching_xref,'A result was pulled from the DB');
is($matching_xref->accession, $xref->accession, 'Retrieved xref is the same as that which was stored');

is($matching_xref->source_id, $source->source_id, 'Source IDs also match');

is($matching_xref->source->name, 'RefSeq', 'Foreign "key" relation works');


# Test auto-populate mechanisms for source, source_url and species tables
$db->populate_metadata($Bin . '/test-data/xref_config.ini');

my $count = $db->schema->resultset('Source')->count(
  { name => 'VGNC' }
);

cmp_ok($count , '==', 1,'One VGNC source in place of four in the original DB');
ok( $db->schema->resultset('Source')->find( { name => 'VGNC' } ), 'Source found' );

# We do not use URLs in the current pipeline, so it's difficult and fruitless to test
# the source_url relationship

my @species = $db->schema->resultset('Species')->search();
is_deeply(
  [map {$_->name} @species],
  [qw/ciona_intestinalis vertebrates danio_rerio xenopus_tropicalis
      pan_troglodytes homo_sapiens canis_familiaris equus_caballus
      bos_taurus mus_musculus rattus_norvegicus/],
  'Species names were all stored');


done_testing;

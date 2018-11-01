use strict;
use warnings;

use Test::More;
use Test::Exception;
use Xref::Test::TestDB;


my $db = Xref::Test::TestDB->new();

$db = Xref::Test::TestDB->new(config_file => 'not_here');
throws_ok { $db->config } qr/does not exist!/,'TestDB grumbles differently when an invalid config file is offered';

# This auto-deploys the schema
$db = Xref::Test::TestDB->new(
  config => { 
    driver => 'SQLite',
    file => 'test.db'
  }
);

ok($db, 'TestDB ready to go');

my $source = $db->schema->resultset('Source')->create({
  name => 'RefSeq',
  status => 'KNOWN',
  source_release => '38',
  download => 'Y',
  priority => 1,
  priority_description => 'Like a boss',
  ordered => 10
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

done_testing;
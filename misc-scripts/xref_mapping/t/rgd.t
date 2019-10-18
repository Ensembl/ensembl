=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use warnings;

use Test::More;

use FindBin '$Bin';

use Xref::Test::TestDB;

use_ok('XrefParser::RGDParser');

my $db = Xref::Test::TestDB->new();

# Initialise two sources:
my $rgd_source = $db->create_db_row('Source',{
  name => 'RGD',
  priority_description => 'Test RGD source',
  priority => 10
});

my $refseq_source = $db->create_db_row('Source',{
  name => 'refseq',
  priority_description => 'Test RefSeq source for dependent xref',
  priority => 20
});

# Create a tame RefSeq xref to which the RGD test data refers
my $refseq_xref = $db->create_db_row('Xref',{
  accession => 'NM1234',
  version => 1,
  label => 'RefSeq innit?',
  source_id => $refseq_source->source_id,
  species_id => 34,
  info_type => 'COORDINATE_OVERLAP',
  info_text => '',
});


my $parser = XrefParser::RGDParser->new($db->dbh);

isa_ok($parser, 'XrefParser::RGDParser');

# Test without any prior RefSeq entries
$parser->run({
  source_id => $rgd_source->source_id,
  species_id => 34,
  files => [ "$Bin/test-data/rgd.txt" ],
});


ok ($db->schema->resultset('Xref')->check_direct_xref({
  accession => 1594427,
  label => '2331ex4-5',
  description => 'class I gene fragment 2331',
}), 'Sample rat Xref has been inserted');

cmp_ok($db->schema->resultset('Xref')->count(), '==', 4, 'Three xrefs in source file and one precursor are now in DB');

# Test that an RGD dependent xref has been assigned to $refseq_xref above

my $hit = $db->schema->resultset('DependentXref')->fetch_dependent_xref('NM1234','1500000');

ok($hit,'Found a dependent Xref');
is ($hit->dependent_xref->accession, '1500000', 'RGD xref is dependent on RefSeq accession');

# Test utility functions
my @sorted = $parser->sort_refseq_accessions(qw/XR100 XP400 NP03201 NM1 NM2/);

is_deeply(\@sorted, [qw/NM1 NM2 NP03201 XP400 XR100/], 'RefSeq accessions are sorted into reasonable order of evidence and utility');

@sorted = $parser->sort_refseq_accessions(qw/a b c d e/);
is_deeply(\@sorted, [], 'Unsanctioned accessions are not included in the sorted output');

my $count = $parser->process_synonyms($refseq_xref->xref_id, 'resistance;is;futile');

cmp_ok($count, '==', 3, 'Three synonyms added from one string');
ok( $db->schema->resultset('Synonym')->check_synonym($refseq_xref->xref_id, 'futile'), 'Third synonym was attached to RefSeq Xref');


done_testing();

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
use Test::Warnings qw( allow_warnings );
use Test::Exception;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts



# Get a DBAdaptor to from the test system
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
ok($multi);
my $db = $multi->get_DBAdaptor("core");
ok($db);

# Should get meaningful type back
debug("get biotype adaptor");
my $biotype_adaptor = $db->get_BiotypeAdaptor();
ok($biotype_adaptor->isa("Bio::EnsEMBL::DBSQL::BiotypeAdaptor"));

# fetch a protein_coding gene object
debug("fetch gene");
my $ga = $db->get_GeneAdaptor();
my $gene = $ga->fetch_by_stable_id("ENSG00000171456");
ok($gene);

# test biotype related methods on this gene object
debug("gene biotype");
is($gene->{'biotype'}, 'protein_coding');
is($gene->biotype, 'protein_coding');

# test gene biotype object
my $biotype = $gene->biotype;
ok($biotype->isa("Bio::EnsEMBL::Biotype"));
is($biotype->object_type, 'gene', 'Biotype is from Gene object');
is($biotype->name, 'protein_coding', 'Biotype name is protein_coding');
is($biotype->biotype_group, 'coding', 'Biotype group is coding');
is($biotype->so_acc, 'SO:0001217', 'Biotype protein_coding refers to SO:0001217');
throws_ok { $biotype->so_acc('test') } qr/so_acc must be a Sequence Ontology accession/,  'so_acc() requires a SO acc like string';
throws_ok { $biotype->object_type('test') } qr/object_type must be gene or transcript/,  'object_type() must be gene or transcript';

# test transcript biotype object
my $transcript = $gene->canonical_transcript;
debug("transcript biotype");
is($transcript->{'biotype'}, 'protein_coding');
is($transcript->biotype, 'protein_coding');
$biotype = $transcript->biotype;
ok($biotype->isa("Bio::EnsEMBL::Biotype"));
is($biotype->object_type, 'transcript', 'Biotype is from Transcript object');
is($biotype->name, 'protein_coding', 'Biotype name is protein_coding');
is($biotype->biotype_group, 'coding', 'Biotype group is coding');
is($biotype->so_acc, 'SO:0000234', 'Biotype protein_coding refers to SO:0000234');

# test fetch biotypes of object_type gene
my $biotypes = $biotype_adaptor->fetch_all_by_object_type('gene');
is(ref $biotypes, 'ARRAY', 'got an array');
is(scalar @{$biotypes}, '2', 'of size 2');

done_testing();

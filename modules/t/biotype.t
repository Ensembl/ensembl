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
use Test::Warnings qw( warning );
use Test::Exception;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts


# Get a DBAdaptor to from the test system
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
ok($multi, "Test DB loaded successfully");
my $db = $multi->get_DBAdaptor("core");
ok($db, "Core DB adaptor loaded successfully");

# Should get meaningful type back
debug("get biotype adaptor");
my $biotype_adaptor = $db->get_BiotypeAdaptor();
ok($biotype_adaptor->isa("Bio::EnsEMBL::DBSQL::BiotypeAdaptor"), "Biotype adaptor loaded successfully");

# fetch a protein_coding gene object
debug("fetch gene");
my $ga = $db->get_GeneAdaptor();
my $gene = $ga->fetch_by_stable_id("ENSG00000171456");
ok($gene, "Gene object loaded successfully");

# test gene biotype object
debug("gene biotype");
is($gene->biotype, 'protein_coding', "Gene biotype is protein_coding");
my $biotype1 = $gene->get_Biotype;
ok($biotype1->isa("Bio::EnsEMBL::Biotype"), "Biotype object retrieved successfully");
is($biotype1->object_type, 'gene', 'Biotype is from Gene object');
is($biotype1->name, 'protein_coding', 'Biotype name is protein_coding');
is($biotype1->biotype_group, 'coding', 'Biotype group is coding');
is($biotype1->so_acc, 'SO:0001217', 'Biotype protein_coding refers to SO:0001217');
throws_ok { $biotype1->so_acc('test') } qr/so_acc must be a Sequence Ontology accession/, 'so_acc() requires a SO acc like string';
throws_ok { $biotype1->object_type('test') } qr/object_type must be gene or transcript/, 'object_type() must be gene or transcript';

# test transcript biotype object
my $transcript = $gene->canonical_transcript;
debug("transcript biotype");
is($transcript->biotype, 'protein_coding', "Trancript biotype is protein_coding");
my $biotype2 = $transcript->get_Biotype;
ok($biotype2->isa("Bio::EnsEMBL::Biotype"), "Biotype object retrieved successfully");
is($biotype2->object_type, 'transcript', 'Biotype is from Transcript object');
is($biotype2->name, 'protein_coding', 'Biotype name is protein_coding');
is($biotype2->biotype_group, 'coding', 'Biotype group is coding');
is($biotype2->so_acc, 'SO:0000234', 'Biotype protein_coding refers to SO:0000234');
ok($transcript->set_Biotype('new_biotype'), "Can successfully set new_biotype");
throws_ok { $gene->set_Biotype() } qr/No argument provided/, 'set_Biotype() requires an argument';

# set biotype with database term
debug("set biotype with db term");
ok($gene->set_Biotype('tRNA'), "Can successfully set biotype to tRNA");
my $biotype3 = $gene->get_Biotype;
ok($biotype3->isa("Bio::EnsEMBL::Biotype"), "Biotype object retrieved successfully");
is($biotype3->object_type, 'gene', 'Biotype is from Gene object');
is($biotype3->name, 'tRNA', 'Biotype name is tRNA');
is($biotype3->biotype_group, 'snoncoding', 'Biotype group is snoncoding');
is($biotype3->so_acc, 'SO:0001263', 'Biotype tRNA refers to SO:0001263');

# set biotype with term not in database
debug("set biotype with term not in db");
ok($gene->set_Biotype('dummy'), "Can successfully set biotype to dummy");
my $biotype4 = $gene->get_Biotype;
ok($biotype4->isa("Bio::EnsEMBL::Biotype"), "Biotype object retrieved successfully");
is($biotype4->object_type, 'gene', 'Biotype is from Gene object');
is($biotype4->name, 'dummy', 'Biotype name is dummy');
is($biotype4->biotype_group, undef, 'Biotype group is not set');
is($biotype4->so_acc, undef, 'Biotype SO acc is not set');
throws_ok { $gene->set_Biotype() } qr/No argument provided/, 'set_Biotype() requires an argument';

# test fetch biotypes of object_type gene
debug("fetch biotypes by object_type");
my $biotypes1 = $biotype_adaptor->fetch_all_by_object_type('gene');
is(ref $biotypes1, 'ARRAY', 'Got an array');
is(scalar @{$biotypes1}, '2', 'of size 2');
is_deeply($biotypes1, [$biotype1, $biotype3], 'with the correct objects');
my $warning1 = warning { 
  $biotypes1 = $biotype_adaptor->fetch_all_by_object_type('none') };
like( $warning1,
    qr/No objects retrieved. Check if object_type 'none' is correct./,
    "Got a warning from fetch_all_by_object_type",
) or diag 'Got warning: ', explain($warning1);
is(ref $biotypes1, 'ARRAY', 'Got an array');
is(scalar @{$biotypes1}, '0', 'of size 0');
is_deeply($biotypes1, [], 'totally empty');
my $biotypes2 = $biotype_adaptor->fetch_all_by_group_object_db_type('coding', 'gene');
is_deeply($biotypes2, [$biotype1], 'fetch_all_by_group_object_db_type retrieves correct data');
my $warning2 = warning { $biotypes2 = $biotype_adaptor->fetch_all_by_group_object_db_type('coding', 'none', 'core') };
like( $warning2,
    qr/No objects retrieved. Check if object_type 'none' is correct./,
    "Got a warning from fetch_all_by_group_object_db_type",
) or diag 'Got warning: ', explain($warning2);

# Test fetching by biotype name
debug("fetch biotypes by name");
my $biotypes3 = $biotype_adaptor->fetch_all_by_name('protein_coding');
is(ref $biotypes3, 'ARRAY', 'Got an array');
cmp_ok(scalar(@$biotypes3), ">", "0", 'Array is not empty');
is_deeply($biotypes3, [$biotype1, $biotype2], 'with the correct objects');
my $warning3 = warning {
  $biotypes3 = $biotype_adaptor->fetch_all_by_name('none') };
like( $warning3,
    qr/No objects retrieved. Check if name 'none' is correct./,
    "Got a warning from fetch_all_by_name",
) or diag 'Got warning: ', explain($warning3);
is(ref $biotypes3, 'ARRAY', 'Got an array');
is(scalar @{$biotypes3}, '0', 'of size 0');
is_deeply($biotypes3, [], 'totally empty');

my $biotypes4 = $biotype_adaptor->fetch_all_by_name('protein_coding', 'gene');
my $biotypes5 = $biotype_adaptor->fetch_by_name_object_type('protein_coding', 'gene');
is_deeply($biotypes4->[0], $biotypes5, "fetch_all_by_name with an object type is identical to fetch_by_name_object_type");

done_testing();

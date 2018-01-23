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

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DnaDnaAlignFeature;

# switch on the note prints
our $verbose = 0;

note("Startup test");
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('ontology');

my $odb = $multi->get_DBAdaptor("ontology");
note("Ontology database instatiated");
ok($odb);

my $human = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $human->get_DBAdaptor("core");
note("Test database instatiated");
ok($db);
my $go_adaptor = $odb->get_OntologyTermAdaptor();

my $accession = "GO:0000217";
my $term = $go_adaptor->fetch_by_accession($accession);
is($term, undef, "GO:0000217 does not exist in the non obsolete list");
$term = $go_adaptor->fetch_by_accession($accession, 1);
ok($term->is_obsolete, "GO:0000217 is obsolete");

$accession = "GO:0003698";
$term = $go_adaptor->fetch_by_accession($accession);
ok($term->name, "GO:0003698 alt_id was fetched");

ok($term->ontology_version() eq 'releases/2016-03-30', "Found ontology version by accession");

$accession = "GO:0003677";
$term = $go_adaptor->fetch_by_accession($accession, 1);
ok(!$term->is_obsolete, "GO:0003677 is not obsolete");
ok(!$term->is_root, "GO:0003677 is not a root");

my $gene;
my $ga = $db->get_GeneAdaptor();

my $genes = $ga->fetch_all_by_GOTerm($term);
is(@{$genes}, 2, "Genes match the GO term");

my $pattern = '%binding%';
my $terms = $go_adaptor->fetch_all_by_name($pattern);
is(@{$terms}, 134, "Found binding terms");

$terms = $go_adaptor->fetch_all_by_name($pattern, undef, 1);
is(@{$terms}, 138, "Found binding terms, including obsolete ones");

ok($term->ontology_version() eq 'releases/2016-03-30', "Found ontology version by name");

my $roots = $go_adaptor->fetch_all_roots();
is(@{$roots}, 2, "Found roots");

my $go_roots = $go_adaptor->fetch_all_roots('GO');
is(@{$go_roots}, 1, "Found go roots");

my $efo_roots = $go_adaptor->fetch_all_roots('EFO');
is(@{$efo_roots}, 0, "Found no efo roots");

#Now go back to the OntologyXref & see if we can do the reverse lookup
{
  my $go_db_links = $genes->[0]->get_all_DBLinks('GO');
  foreach my $dbentry (@{$go_db_links}) {
    my $term = $dbentry->get_OntologyTerm();
    my $direct_term = $go_adaptor->fetch_by_accession($dbentry->primary_id());
    is_deeply($term, $direct_term, 'Fetching the OntologyTerm from the OntologyXref should match the same object from OntologyTermAdaptor');
  }
}
$term = $go_adaptor->fetch_by_accession('GO:0000182',1); # unintentionally picked an obsolete term for testing on.
my $term_list = $go_adaptor->fetch_all_by_descendant_term($term);
my $inclusive_term_list = $go_adaptor->fetch_all_by_descendant_term($term,undef,undef,1);
my $other_term_list = $term->ancestors();
ok (scalar(@$term_list) == scalar(@$inclusive_term_list) - 1, "Zero_distance flag on fetch_all_by_descendant_term");
is(scalar(@$term_list), scalar(@$other_term_list), "Fetching all by descendant is the same as fetching all ancestors for term");
is($term_list->[0]->accession, $other_term_list->[0]->accession, "Fetching all by descendant is the same as fetching all ancestors for term");

my $parent_list = $go_adaptor->fetch_all_by_parent_term($term);
my $other_parent_list = $term->children();
is(scalar(@$parent_list), 2, 'Term has 2 parent');
is(scalar(@$parent_list), scalar(@$other_parent_list), "Fetching all by parent is the same as fetching all children for term");
is($parent_list->[0], $other_parent_list->[0], "Same terms returned");

my $child_list = $go_adaptor->fetch_all_by_child_term($term);
my $other_child_list = $term->parents();
is(scalar(@$child_list), 1, 'Term has 1 child');
is(scalar(@$child_list), scalar(@$other_child_list), "Fetching all terms by child is the same as fetching all parents for term");
is($child_list->[0], $other_child_list->[0], "Fetching all terms by child is the same as fetching all parents for term");

my $chart = $go_adaptor->_fetch_ancestor_chart($term, 'GO');
ok(%$chart, 'Can fetch ancestor chart');

my $new_term = $go_adaptor->fetch_by_dbID($term->dbID, 1);
is($new_term->accession, $term->accession, "Fetched the same term using dbID");

my @dbid_list = ($term->dbID, $term_list->[0]->dbID);
my $list = $go_adaptor->fetch_all_by_dbID_list(\@dbid_list);
is($list->[0]->accession, $term_list->[0]->accession, "Fetched the correct term using list of dbIDs");
my $obsolete_list = $go_adaptor->fetch_all_by_dbID_list(\@dbid_list, 1);
is($obsolete_list->[0]->accession, $term->accession, "First dbID is obsolete");
is(scalar(@$list), scalar(@$obsolete_list) - 1, "One obsolete term found");

my $alts = $go_adaptor->fetch_all_alt_ids('GO:0000182');
is(scalar(@$alts), 0, "No alternative accessions for GO:0000182");

my $all = $go_adaptor->fetch_all();
my $all_obsolete = $go_adaptor->fetch_all(1);
is(scalar(@$all), 164, "164 terms found");
is(scalar(@$all_obsolete), 168, "168 terms found when including obsolete ones");

done_testing();

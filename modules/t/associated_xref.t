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
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;

note( "Startup test" );
#
# 1 Test started
#
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );

note( "Test database instantiated" );
#
# 2 Database instatiated
#
ok( $db );

my $dbEntryAdaptor = $db->get_DBEntryAdaptor();

#
# 3 Get a gene to attach the term to
#
my $ga  = $db->get_adaptor('Gene');
my $gene = $ga->fetch_by_stable_id("ENSG00000171456");
debug("Gene->fetch_by_stable_id()");
ok($gene);

$multi->hide("core", 'xref', 'object_xref', 'associated_xref', 'associated_group');

#
# 4 Loading ontology terms and adding associate_xrefs
#
my $ont_xref = Bio::EnsEMBL::OntologyXref->new
  (
   -primary_id => "GO:0000001", 
   -dbname => "GO",
   -release => "1", 
   -display_id => "Test GO term",
   -description => "Amuch longer description of the GO term",
   -info_type => "DIRECT",
   -analysis => undef,
   );   

my $pmid_xref = Bio::EnsEMBL::DBEntry->new (
                  -primary_id => "123456789",
                  -dbname     => "PUBMED",
                  -release    => "1",
                  -display_id => "PMID:123456789",
                  -description => "A paper",
                  -info_type  => "DIRECT"
                );

my $during_xref = Bio::EnsEMBL::DBEntry->new (
                  -primary_id => "GO:0007067",
                  -dbname     => "GO",
                  -release    => "1",
                  -display_id => "mitotic nuclear division",
                  -description => "A cell cycle process comprising the steps by which the nucleus of a eukaryotic cell divides; the process involves condensation of chromosomal DNA into a highly compacted form. Canonically, mitosis produces two daughter nuclei whose chromosome complement is identical to that of the mother cell.",
                  -info_type  => "DIRECT"
                );

my $with_xref = Bio::EnsEMBL::DBEntry->new (
                  -primary_id => "Q9BR19",
                  -dbname     => "Uniprot/SPTREMBL",
                  -release    => "1",
                  -display_id => "Q9BR19",
                  -info_type  => "DIRECT"
                );

$ont_xref->add_linkage_type( "IC", $pmid_xref ); # Linkage type on own

$ont_xref->add_associated_xref(
  $during_xref,
  $pmid_xref,
  "during",
  1,
  1
);

$ont_xref->add_associated_xref(
  $with_xref,
  $pmid_xref,
  "with",
  1,
  2
);


my $ont_xref_id = $dbEntryAdaptor->store($ont_xref, $gene->dbID, "Gene");
note("Xref_id from insert: ".$ont_xref_id);


#
# HC1 - Check that there are 2 entries in the associated_xref table
#
my $assoc_count = $db->dbc->prepare( 'select count(*) from associated_xref' );
#$assoc_count->bind_param(1, $ont_xref_id);
$assoc_count->execute();

my ( $assoc_xref_count )  = $assoc_count->fetchrow_array();
$assoc_count->finish();
is($assoc_xref_count, 2, "Count in the associated_xref table.");


#
# HC2 - Check that there is a single entry in the associated_group table
#
$assoc_count = $db->dbc->prepare( 'select count(*) from associated_group' );
$assoc_count->execute();

( $assoc_xref_count )  = $assoc_count->fetchrow_array();

is($assoc_xref_count, 1, "Count in the associated_group table.");
$assoc_count->finish();


#
# 5 Retrieve associated_xref entries from the db.
#
my $gene2 = $ga->fetch_by_stable_id("ENSG00000171456");
debug("Gene->fetch_by_stable_id()");
ok($gene2, "Found the required gene");

my @gene_DBentries = @{ $gene2->get_all_DBEntries() };
my $ont_xref_retrieved = 0;
my $ont_xref_hasAssoc  = 0;
foreach my $dbentry (@gene_DBentries) {
  if (ref $dbentry eq 'Bio::EnsEMBL::OntologyXref') {
    if ($dbentry->dbname eq 'GO' and $dbentry->primary_id eq 'GO:0000001') {
      $ont_xref_retrieved = 1;
      my $annot_ext = $dbentry->get_all_associated_xrefs();
      $ont_xref_hasAssoc  = scalar(keys %{ $annot_ext });
    }
  }
}

#
# HC3 - Test the retrieval of the ontology term with paired asociated_xrefs
#
is($ont_xref_retrieved, 1, "Ontology DBEntry has been retrieved");
is($ont_xref_hasAssoc, 1, "Count of the number associated xrefs");

$multi->restore;

done_testing();

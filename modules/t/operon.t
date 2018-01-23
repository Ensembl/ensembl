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
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Operon;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::DBSQL::OperonAdaptor;
use Bio::EnsEMBL::DBSQL::OperonTranscriptAdaptor;
debug("Startup test");
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $dba = $multi->get_DBAdaptor("core");

debug("Test database instatiated");
ok($dba);

# get a slice
my $slice = $dba->get_SliceAdaptor()->fetch_by_seq_region_id(469283);

# create operon
my $start         = 31225346;
my $end           = 31225946;
my $strand        = 1;
my $display_label = "accBC";
my $analysis = $dba->get_AnalysisAdaptor->fetch_by_logic_name("Genscan");
ok(defined $analysis);
my $operon = Bio::EnsEMBL::Operon->new(
	-START         => $start,
	-END           => $end,
	-STRAND        => $strand,
	-SLICE         => $slice,
	-DISPLAY_LABEL => $display_label,
	-ANALYSIS => $analysis);
$operon->add_DBEntry( Bio::EnsEMBL::DBEntry->new( -DBNAME     => 'EMBL',
												  -RELEASE    => 1,
												  -PRIMARY_ID => 'XZY',
												  -DISPLAY_ID => '123',
												  -ANALYSIS   => $analysis ) );

# check it has the correct properties
is( $display_label, $operon->display_label(),     "Operon name" );
is( $start,         $operon->seq_region_start(),  "Operon start" );
is( $end,           $operon->seq_region_end(),    "Operon end" );
is( $strand,        $operon->seq_region_strand(), "Operon strand" );
is( $analysis,        $operon->analysis(), "Analysis" );

my $operon_adaptor = Bio::EnsEMBL::DBSQL::OperonAdaptor->new($dba);

# store operon
$operon_adaptor->store($operon);
ok( defined $operon->dbID() );
# retrieve operon
my $operon2 = $operon_adaptor->fetch_by_dbID( $operon->dbID() );
is( $operon2->dbID(),             $operon->dbID(),             "Operon ID" );
is( $operon2->display_label(),    $operon->display_label(),    "Operon name" );
is( $operon2->seq_region_start(), $operon->seq_region_start(), "Operon start" );
is( $operon2->seq_region_end(),   $operon->seq_region_end(),   "Operon end" );
is( $operon2->seq_region_strand(),
	$operon->seq_region_strand(),
	"Operon strand" );
is( $operon2->analysis(),
	$operon->analysis(),
	"Analysis" );

SKIP: {
  skip 'No registry support for SQLite yet', 1 if $operon_adaptor->dbc->driver() eq 'SQLite';

  #test the get_species_and_object_type method from the Registry
  my $registry = 'Bio::EnsEMBL::Registry';
  my ( $species, $object_type, $db_type ) = $registry->get_species_and_object_type('16152-16153-4840');
  ok( $species eq 'homo_sapiens' && $object_type eq 'Operon');
}

debug ("Operon->list_stable_ids");
my $stable_ids = $operon_adaptor->list_stable_ids();
ok (@{$stable_ids});


$operon = $operon_adaptor->fetch_by_stable_id('16152-16153-4840');
debug( "Operon->fetch_by_stable_id()" );
ok( $operon );

$operon->stable_id_version('16152-16153-4841.4');
is($operon->stable_id, '16152-16153-4841', 'Stable id set with stable_id_version');
is($operon->version, 4, 'Version set with stable_id_version');
is($operon->stable_id_version, '16152-16153-4841.4', 'Stable id and version from stable_id_version');

$operon->stable_id_version('16152-16153-4842');
is($operon->stable_id, '16152-16153-4842', 'Stable id set with stable_id_version');
is($operon->version, undef, 'Version undef from stable_id_version');
is($operon->stable_id_version, '16152-16153-4842', 'Stable id and no version from stable_id_version');

$operon = $operon_adaptor->fetch_by_stable_id('16152-16153-4840.1');
ok($operon->stable_id eq '16152-16153-4840', 'fetch_by_stable_id with version');

$operon = $operon_adaptor->fetch_by_stable_id('16152-16153-4840.1a');
ok(! defined($operon), 'fetch_by_stable_id with bad version');

$operon = $operon_adaptor->fetch_by_stable_id_version('16152-16153-4840', 1);
ok($operon->stable_id eq '16152-16153-4840', 'fetch_by_stable_id_version with version');

$operon = $operon_adaptor->fetch_by_stable_id_version('16152-16153-4840', '1a');
ok(! defined($operon), 'fetch_by_stable_id_version with bad version');

#19
my @operons = @{ $operon_adaptor->fetch_all_versions_by_stable_id('16152-16153-4840') };
debug("fetch_all_versions_by_stable_id");
ok( scalar(@operons) == 1 );


#20

$operon = $operon_adaptor->fetch_by_operon_transcript_stable_id('T16152-16153-4840');
debug( "Operon->fetch_by_operon_transcript_stable_id()" );
ok( $operon );


#21-24

#
# Operon remove test
#

$multi->save( "core", "operon", "operon_transcript",
	      'xref', "object_xref", "ontology_xref", "identity_xref", 'meta_coord');

$operon = $operon_adaptor->fetch_by_stable_id( "16152-16153-4840" );

my $operon_count = count_rows( $dba, "operon" );
my $operon_trans_count = count_rows( $dba, "operon_transcript" );

my $ots = scalar( @{$operon->get_all_OperonTranscripts() } );

debug( "Operons before ".$operon_count );
debug( "OperonTranscripts before ".$operon_trans_count );
debug( "Operon has ".$ots." transcripts" );

$operon_adaptor->remove( $operon );

ok( count_rows( $dba, "operon" ) == ( $operon_count - 1 ));
ok( count_rows( $dba, "operon_transcript" ) == ($operon_trans_count-$ots));

ok(!defined($operon->dbID()));
ok(!defined($operon->adaptor()));

$multi->restore('core');

done_testing();

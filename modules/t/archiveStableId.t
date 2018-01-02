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
no warnings qw(uninitialized);

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0;


#
# 1 ArchiveStableId adaptor compiles
#
ok(1);


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db    = $multi->get_DBAdaptor('core');

my $asia = $db->get_ArchiveStableIdAdaptor();


#
# 2-4 ArchiveStableId retrieval
#
my $asi = $asia->fetch_by_stable_id("T1");
is( $asi->release, 2, "T1 is from release 2");

$asi = $asia->fetch_by_stable_id_version("T2", 3);
is( $asi->release, 3, "T2 is from release 3");

$asi = $asia->fetch_by_stable_id_dbname("T1", "release_2");
is( $asi->release, 2, "T1 is from release 2");


#
# 5 retrieval of an archiveStableId
#
$asi = $asia->fetch_by_stable_id( "G1" );
_print_asi( $asi );

ok( $asi );


#
# 6 retrieval of the event related to a specific stable id
#
my $event = $asi->get_event("G2");

is(ref($event), 'Bio::EnsEMBL::StableIdEvent', "A stable id event was fetched");
is(sprintf("%.6f", $event->score), sprintf("%.6f", 0.54), "Mapping score between G1 and G2");
my $string = $event->ident_string();

my $old_archive_stable_id = $event->old_ArchiveStableId;
my $new_archive_stable_id = $event->new_ArchiveStableId;

like($string, qr/G2.3 \(3\) -> G1.2 \(4\) \[0.54/, "Event string");

is($new_archive_stable_id, $asi, "Initial archive is new archive");
is($old_archive_stable_id->stable_id, "G2", "Old stable id");
is($new_archive_stable_id->stable_id, "G1", "New stable id");

$event = $old_archive_stable_id->get_event("G1");

is(sprintf("%.6f", $event->score), sprintf("%.6f", 0.54), "Mapping score between G1 and G2");

$old_archive_stable_id = $event->old_ArchiveStableId;
$new_archive_stable_id = $event->new_ArchiveStableId;

is($old_archive_stable_id->stable_id, "G2", "Old stable id");
is($new_archive_stable_id->stable_id, "G1", "New stable id");


#
# 7 how many predecessors does it have
#
my $pre_asis = $asi->get_all_predecessors();
is( scalar( @$pre_asis ), 2, "G1 has 2 predecessors" );

for my $asi ( @$pre_asis ) {
  debug( "\tPre G1" );
  _print_asi( $asi );
}


#
# 8 transcripts for a gene
#
my $transcripts = $pre_asis->[0]->get_all_transcript_archive_ids();

for my $asi ( @$transcripts ) {
  debug( "\tTranscripts G1" );
  _print_asi( $asi );
  
  my $tl = $asi->get_all_translation_archive_ids();
  foreach my $asi2 (@$tl) {
    _print_asi( $asi2 );
  }
}

is( scalar( @$transcripts ), 1, "G1 has 1 transcript");


#
# 9 no predecessor case
#
$pre_asis = $pre_asis->[0]->get_all_predecessors();
debug( "\tPredecessors: ".scalar( @$pre_asis ) );

is( scalar( @$pre_asis ), 0, "No predecessors found" );


#
# 10 successor case
#
$asi = $asia->fetch_by_stable_id_dbname( "G4", "release_1" );
my $succ_asis = $asi->get_all_successors();
 
for my $asi ( @$succ_asis ) {
  debug( "\tSucc G4.1" );
  _print_asi( $asi );
}

is( scalar( @$succ_asis ), 1, "G4 has 1 sucessor" );

#
# 11 no successor case
#
$succ_asis = $succ_asis->[0]->get_all_successors();

for my $asi ( @$succ_asis ) {
  debug( "\tSucc Succ G4.1" );
  _print_asi( $asi );
}

is( scalar( @$succ_asis ), 0, "G4.1 has no sucessors");


#
# 12 fetch_successor_history
#
$asi = $asia->fetch_by_stable_id_dbname( "G2", "release_1" );
my $asis = $asia->fetch_successor_history( $asi );

debug( "\tCurrently related from G2.release_1" );
for my $asi ( @$asis ) {
 _print_asi( $asi );
}

is( $asis->[-1]->db_name, "release_4", "Current release for G2 is release 4");
is( scalar @$asis, 5, "G2 has 5 sucessors");

#
# 13-17 history tree
#
$asi = $asia->fetch_by_stable_id_dbname( "G2", "release_1" );
my $history = $asi->get_history_tree;

my @asis = @{ $history->get_all_ArchiveStableIds };
is( scalar(@asis), 9, "G2 history has 9 related archives");

my @events = @{ $history->get_all_StableIdEvents };
is( scalar(@events), 10, "G2 history has 10 related events");

is( scalar(@{ $history->get_release_display_names }), 4, "G2 has 4 display names");
is( scalar(@{ $history->get_unique_stable_ids }), 3, "G2 has 3 unique stable ids");

my ($x, $y) = @{ $history->coords_by_ArchiveStableId($asi) };
ok( $x == 0 and $y == 1 );


#
# 18-19 check for current version and fetch latest incarnation
#
ok( ! $asi->is_latest, 'Not on the latest version so is_latest is false');

$asi = $asi->get_latest_incarnation;
ok($asi->is_latest(), 'Latest incarnation must be the latest version');
is($asi->version, 4, 'Latest version is 4');

#
# 20 associated IDs in archive
#
$asi = $asia->fetch_by_stable_id_version( "G2", "2" );
my @assoc = @{ $asi->get_all_associated_archived };
ok( scalar(@assoc) == 2 and
    $assoc[0]->[0]->type eq 'Gene' and
    $assoc[0]->[1]->type eq 'Transcript' and
    $assoc[0]->[2]->type eq 'Translation' and
    $assoc[0]->[3] =~ /^PT/
);


#
# 21 archived peptide sequence
#
$asi = $asia->fetch_by_stable_id_version("P2", 1);
ok( $asi->get_peptide eq 'PTWOVERSIONONE*' );

#
# Test looking up active ids
#
my $archive_obj = $asia->fetch_by_stable_id('ENSG00000171456');
is($archive_obj->stable_id, 'ENSG00000171456', 'fetch_by_stable_id with active stable_id');
ok($archive_obj->is_current, 'Is the current stable_id');

$archive_obj = $asia->fetch_by_stable_id('ENSG00000171456.1');
is($archive_obj->stable_id, 'ENSG00000171456', 'fetch_by_stable_id with active stable_id with version');
ok($archive_obj->is_latest, 'Is the latest stable_id');

#
# debug helper
#
sub _print_asi {
  my $asi = shift;

  debug( "\ttype: ".$asi->type().
         "\n\tstable id: ".$asi->stable_id().
	 "\n\tversion: ".$asi->version().
	 "\n\tdbname: ".$asi->db_name().
	 "\n\tTranscripts: ".(join(", ", map { $_->stable_id } @{ $asi->get_all_transcript_archive_ids })).
	 "\n\tTranslations: ".(join(", ", map { $_->stable_id } @{ $asi->get_all_translation_archive_ids })).
	 "\n\tPeptide: ".$asi->get_peptide."\n" );
}
 
done_testing();

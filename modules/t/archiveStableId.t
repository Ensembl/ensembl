use strict;
use warnings;
no warnings qw(uninitialized);

BEGIN { $| = 1;  
	use Test;
	plan tests => 20;
}

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
ok( $asi->release == 2);

$asi = $asia->fetch_by_stable_id_version("T2", 3);
ok( $asi->release == 3);

$asi = $asia->fetch_by_stable_id_dbname("T1", "release_2");
ok( $asi->release == 2);


#
# 5 retrieval of an archiveStableId
#
$asi = $asia->fetch_by_stable_id( "G1" );
_print_asi( $asi );

ok( $asi );


#
# 6 how many predecessors does it have
#
my $pre_asis = $asi->get_all_predecessors();
ok( scalar( @$pre_asis ) == 2 );

for my $asi ( @$pre_asis ) {
  debug( "\tPre G1" );
  _print_asi( $asi );
}


#
# 7 transcripts for a gene
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

ok( scalar( @$transcripts ) == 1);


#
# 8 no predecessor case
#
$pre_asis = $pre_asis->[0]->get_all_predecessors();
debug( "\tPredecessors: ".scalar( @$pre_asis ) );

ok( scalar( @$pre_asis ) == 0 );


#
# 9 successor case
#
$asi = $asia->fetch_by_stable_id_dbname( "G4", "release_1" );
my $succ_asis = $asi->get_all_successors();
 
for my $asi ( @$succ_asis ) {
  debug( "\tSucc G4.1" );
  _print_asi( $asi );
}

ok( scalar( @$succ_asis ) == 1 );

#
# 10 no successor case
#
$succ_asis = $succ_asis->[0]->get_all_successors();

for my $asi ( @$succ_asis ) {
  debug( "\tSucc Succ G4.1" );
  _print_asi( $asi );
}

ok( scalar( @$succ_asis ) == 0 );


#
# 11 fetch_successor_history
#
$asi = $asia->fetch_by_stable_id_dbname( "G2", "release_1" );
my $asis = $asia->fetch_successor_history( $asi );

debug( "\tCurrently related from G2.release_1" );
for my $asi ( @$asis ) {
 _print_asi( $asi );
}

ok(( $asis->[-1]->db_name eq "release_4" ) &&
   ( scalar @$asis == 5 ));

#
# 12-16 history tree
#
$asi = $asia->fetch_by_stable_id_dbname( "G2", "release_1" );
my $history = $asi->get_history_tree;

my @asis = @{ $history->get_all_ArchiveStableIds };
ok( scalar(@asis) == 9);

my @events = @{ $history->get_all_StableIdEvents };
ok( scalar(@events) == 10);

ok( scalar(@{ $history->get_release_display_names }) == 4);
ok( scalar(@{ $history->get_unique_stable_ids }) == 3);

my ($x, $y) = @{ $history->coords_by_ArchiveStableId($asi) };
ok( $x == 0 and $y == 1 );


#
# 17-18 check for current version and fetch latest incarnation
#
ok( ! $asi->is_latest );

$asi = $asi->get_latest_incarnation;
ok( $asi->is_latest and $asi->version == 3 );

#
# 19 associated IDs in archive
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
# 20 archived peptide sequence
#
$asi = $asia->fetch_by_stable_id_version("P2", 1);
ok( $asi->get_peptide eq 'PTWOVERSIONONE*' );


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
  

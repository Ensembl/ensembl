use strict;
use warnings;
no warnings qw(uninitialized);

BEGIN { $| = 1;  
	use Test;
	plan tests => 9;
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

my $asi = $asia->fetch_by_stable_id( "G1" );

_print_asi( $asi );

#
# 2 retrieval of an archiveStableId
#
ok( $asi );

my $pre_asis = $asi->get_all_predecessors();

for my $asi ( @$pre_asis ) {
  debug( "\tPre G1" );
  _print_asi( $asi );
}


#
# 3 how many predecessors does it have
#
ok( scalar( @$pre_asis ) == 2 );

my $transcripts = $pre_asis->[0]->get_all_transcript_archive_ids();

for my $asi ( @$transcripts ) {
  debug( "\tTranscripts G1" );
  _print_asi( $asi );
  
  #get_translation_archive_id was changed to give back listref.
  #this makes the function poorly named, but it is what the
  #webteam uses so....
  my $tl = $asi->get_all_translation_archive_ids();
  foreach my $asi2 (@$tl) {
    _print_asi( $asi2 );
  }
}

#
# 4 transcripts for a gene
#
ok( scalar( @$transcripts ) == 1);



$pre_asis = $pre_asis->[0]->get_all_predecessors();
debug( "\tPredecessors: ".scalar( @$pre_asis ) );


#
# 5 no predecessor case
#
ok( scalar( @$pre_asis ) == 0 );

$asi = $asia->fetch_by_stable_id_dbname( "G4", "release_1" );
my $succ_asis = $asi->get_all_successors();
 
for my $asi ( @$succ_asis ) {
  debug( "\tSucc G4.1" );
  _print_asi( $asi );
}

#
# 6 successor case
#
ok( scalar( @$succ_asis ) == 1 );

$succ_asis = $succ_asis->[0]->get_all_successors();

for my $asi ( @$succ_asis ) {
  debug( "\tSucc Succ G4.1" );
  _print_asi( $asi );
}


#
# 7 no successor case
#

ok( scalar( @$succ_asis ) == 0 );

#
# 8 fetch_successor_history
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
# 9 reject unknown stable ids
#

ok( ! defined $asia->fetch_by_stable_id_dbname( "FooBar", "release_unknown" ));



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
  

use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 34;
}

use MultiTestDB;
use TestUtils qw(debug test_getter_setter);
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::Exon;

our $verbose = 0;

#
# 1 PredictionTranscript compiles
#
ok(1);

my $multi = MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_chr_start_end("20", 30_252_000, 31_252_001 );

my $p_transs = $slice->get_all_PredictionTranscripts();

#
# 2 Verify prediciton transcripts can be obtained
#
ok( scalar( @$p_transs ) );


if($verbose) {
  foreach my $pt (@$p_transs) {
    print $pt->stable_id , "\n";
  }
}

my ($pt) = @$p_transs;

my $exons = $pt->get_all_Exons;

#
# 3 test get all exons
#
ok(scalar @$exons);

#
#  4 test new
#
my $new_pt = new Bio::EnsEMBL::PredictionTranscript(@$exons);
ok(scalar @{$new_pt->get_all_Exons});

#
# 5 test stable_id
#
ok($pt->stable_id =~ /(\w+\.\d+\.\d+\.\d+)\.(\d+)\.(\d+)/);
#my ($ctg, $ctg_start, $ctg_end) = ($1, $2, $3);


#
# 6 test coding start
#
ok(&TestUtils::test_getter_setter($pt, 'coding_region_start', 6));

#
# 7 test coding end
#
ok(&TestUtils::test_getter_setter($pt, 'coding_region_end', 7));


#
# 8 test start
#
ok(&TestUtils::test_getter_setter($pt, 'start', 8));

#
# 9 test end
#
ok(&TestUtils::test_getter_setter($pt, 'end', 9));


#
# 10 test analysis 
#
my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name('Vertrna');
ok(&TestUtils::test_getter_setter($pt, 'analysis', $analysis));

#
# 11 test dbID
#
ok(&TestUtils::test_getter_setter($pt, 'dbID', 11));

#
# 12 test adaptor
#
my $pta = $db->get_PredictionTranscriptAdaptor;
ok(&TestUtils::test_getter_setter($pt, 'adaptor', $pta));

#
# 13-17 test add Exon
#
my $exon = new Bio::EnsEMBL::Exon;
$exon->start($pt->start - 20);
$exon->end($pt->end + 20);

my @old_exons = @{$pt->get_all_Exons};
my $count = scalar(@old_exons);

$pt->add_Exon($exon);
#check that transcript start + end updated
ok($pt->start == $exon->start);
ok($pt->end == $exon->end);
$exons = $pt->get_all_Exons( 1 );
#check that there is one more exon
ok($count + 1 == scalar(@$exons)); 
#check that the last exon is the exon added
ok($exons->[$#$exons] eq $exon); 

$pt->add_Exon($exon, 3);
#check that third exon is exon added
ok($exons->[2] eq $exon); 


#
# 18-22 test flush exons
#
$pt->flush_Exons;
ok(scalar @{$pt->get_all_Exons} == 0);
ok(!defined $pt->start);
ok(!defined $pt->end);
ok(!defined $pt->coding_start);
ok(!defined $pt->coding_end);

#restore old exons
my $pos = 0;
while(@old_exons) {
  $pos++;
  my $e = shift @old_exons;
  $pt->add_Exon($e, $pos);
}

#
# 23 test get_all_translateable_Exons
#
ok(scalar @{$pt->get_all_translateable_Exons} == scalar @{$pt->get_all_Exons});

#
# 24 test sort executes
#
ok($pt->sort || 1);

#
# 25 test get_exon_count
#
ok($pt->get_exon_count == scalar @{$pt->get_all_Exons});

#
# 26 test set_exon_count
#
$count = $pt->get_exon_count;
$pt->set_exon_count(26);
ok(scalar @{$pt->get_all_Exons( 1 )} == 26); # test internal exon array expansion
$pt->set_exon_count($count);

#
# 27 test length
#
my $len = 0;
foreach my $ex (@{$pt->get_all_Exons}) {
  if( defined $ex ) { $len += $ex->length };
}
ok($len == $pt->length);

#
# 28 test translate
#
ok(length $pt->translate->seq);

#
# 29 test get cdna
#
my $cstart = 1;
my $cend   = $pt->length;
ok(length $pt->get_cdna($cstart, $cend));

#
# 30 test pep2genomic
#
my $pend = $pt->length / 3;
my $pstart = 2;

my $defined_exons_count = 0;
foreach my $e (@{$pt->get_all_Exons}) {
  if(defined $e) {
    $defined_exons_count++;
  }
}
# should return genomic coords for each exon since covers entire peptide
ok($defined_exons_count == $pt->pep2genomic($pstart, $pend));


#
# 31 test cdna2genomic
#
ok($defined_exons_count == $pt->cdna2genomic($cstart, $cend));


#
# 32 test type
#
ok(&TestUtils::test_getter_setter($pt, 'type', 'test'));

#
# 33 test fetch_by_stable_id
#

my $stable_id = 'AL031658.11.1.162976.122801.143660';

$pt = $pta->fetch_by_stable_id($stable_id);
ok($pt->stable_id eq $stable_id);

# 34 list_dbIDs
my $ids = $pta->list_dbIDs();
ok (@{$ids});



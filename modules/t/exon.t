use lib 't';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 21;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use MultiTestDB;
use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts

$loaded = 1;
my $multi = MultiTestDB->new();

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

ok($db);


# Exon specific tests

my $exonad = $db->get_ExonAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

my $slice = $slice_adaptor->fetch_by_region('chromosome', '20',
                                            30_811_000,
                                            32_000_000);
ok($exonad);

my $exon = Bio::EnsEMBL::Exon->new();


$exon->start(1000);
ok(&test_getter_setter($exon, 'start', 200));

$exon->end(1400);
ok(&test_getter_setter($exon, 'end', 400));

$exon->strand(1);
ok(&test_getter_setter($exon, 'strand', -1));

$exon->phase(0);
ok(&test_getter_setter($exon, 'phase', -1));

$exon->slice( $slice );
ok(&test_getter_setter($exon, 'slice', $slice));

# should try to store (!)
$exon->end_phase( -1 );
ok(&test_getter_setter($exon, 'end_phase', 1));

#
# find supporting evidence for the exon
#
my @evidence = ();
my @fs = ();
push @fs, @{$db->get_DnaAlignFeatureAdaptor->fetch_all_by_Slice($slice)};
push @fs, @{$db->get_ProteinAlignFeatureAdaptor->fetch_all_by_Slice($slice)};

while(my $f = shift @fs) {
  #debug("feature at: " . $f->start . "-" . $f->end);
  next if $f->start > $exon->end || $f->end < $exon->start;
  push(@evidence, $f);
  # cheat it into storing it again
  $f->dbID( undef );
  $f->adaptor( undef );
}

my $count = scalar(@evidence);
debug("adding $count supporting features");
$exon->add_supporting_features(@evidence);

$multi->hide( "core", "exon", "supporting_feature", 
	      "protein_align_feature", "dna_align_feature");

$exonad->store($exon);

ok($exon->dbID() == 1 && $exon->adaptor == $exonad);

# now test fetch_by_dbID

my $newexon = $exonad->fetch_by_dbID($exon->dbID);

ok($newexon);



debug("exon chr start  = " . $exon->start);
debug("exon chr end    = " . $exon->end);
debug("exon chr strand = " . $exon->strand);

debug("newexon start  = " . $newexon->start());
debug("newexon end    = " . $newexon->end());
debug("newexon strand = " . $newexon->strand());

ok($newexon->start == 30811999 &&
   $newexon->end == 30812399 &&
   $newexon->strand==1);


#
# Test transform to another slice
#
$slice = $exon->slice();
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                         $slice->seq_region_name,
                                         $slice->start + $exon->start - 11,
                                         $slice->start + $exon->end + 9);
$exon = $exon->transfer($slice);
debug("exon start  = " . $exon->start);
debug("exon end    = " . $exon->end);
debug("exon strand = " . $exon->strand);
ok($exon->start == 11 && $exon->end == 411 && $exon->strand==1);


#
# Test Transform to contig coord system
#
$exon = $exon->transform('contig');

debug("exon start  = " . $exon->start);
debug("exon end    = " . $exon->end);
debug("exon strand = " . $exon->strand);
debug("exon seq_region = " . $exon->slice->seq_region_name);

ok($exon->start == 913);
ok($exon->end   == 1313);
ok($exon->strand == 1);
ok($exon->slice->seq_region_name eq 'AL034550.31.1.118873');


#regression test, supporting evidence was lost post transform before...
my $se_count = scalar(@{$exon->get_all_supporting_features});

debug("Got $se_count supporting feature after transform");
ok($se_count == $count);

#make sure that supporting evidencewas stored correctly
$se_count = scalar(@{$newexon->get_all_supporting_features});
debug("Got $se_count from newly stored exon");
ok($se_count == $count);


# list_ functions
debug ("Exon->list_dbIDs");
my $ids = $exonad->list_dbIDs();
ok (@{$ids});

debug ("Exon->list_stable_ids");
my $stable_ids = $exonad->list_stable_ids();
ok (@{$stable_ids});

$multi->restore();


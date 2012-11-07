use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Slice;

our $verbose= 0;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

#
# Test fetch_by_Slice_start_end_strand
#
my $slice_adaptor = $db->get_SliceAdaptor;
my $seq_adaptor = $db->get_CompressedSequenceAdaptor();


my $seq =  'ACTGAAANTTANNNATYTTTAAATTACCC';
my $len = length($seq);

my $contig_cs = $db->get_CoordSystemAdaptor->fetch_by_name('contig');

$multi_db->save('core', 'dnac', 'seq_region');

#we need to create a fake seq region because the length of the seq region
#needs to match the length of the sequence we are inserting. Otherwise
#we can get weird padding at the end due to the last byte being not
#fully packed but still labelled as non-gap 
my $sth = 
  $db->dbc->prepare('INSERT INTO seq_region (name, length, coord_system_id) ' .
               'VALUES (?,?,?)');
$sth->execute('testfrag', $len,$contig_cs->dbID);


my $slice = $slice_adaptor->fetch_by_region('contig', 'testfrag');

my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);


debug("Storing sequence:   $seq");
$seq_adaptor->store($seq_region_id,$seq);

ok(1);

my $new_seq = ${$seq_adaptor->fetch_by_Slice_start_end_strand($slice, 1, 
                                                              $len, 1)};

ok($seq eq $new_seq);


debug("Retrieved sequence: $seq");

my $flanking = 5;
$new_seq = ${$seq_adaptor->fetch_by_Slice_start_end_strand($slice, 
                                                           1-$flanking, 
                                                           $len+$flanking, 1)};


ok(('N' x $flanking) . $seq . ('N' x $flanking) eq $new_seq);

debug("Retrieved sequence (with $flanking flanking): $new_seq");

$multi_db->restore('core', 'dnac', 'seq_region');

done_testing();

use lib 't';
use strict;
use warnings;


BEGIN { $| = 1;  
	use Test;
	plan tests => 18;
}

use MultiTestDB;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use TestUtils qw(test_getter_setter debug);

our $verbose = 0;

my ($CHR, $START, $END, $FLANKING) = ("20", 30_252_000, 31_252_001, 1000);

#
# slice adaptor compiles
#
ok(1);

my $multi = MultiTestDB->new;
my $db    = $multi->get_DBAdaptor('core');


#
# SliceAdaptor::new
#
my $slice_adaptor = Bio::EnsEMBL::DBSQL::SliceAdaptor->new($db->_obj);
ok($slice_adaptor->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor'));
ok($slice_adaptor->db);

#
# fetch_by_region
#
my $slice = $slice_adaptor->fetch_by_region('toplevel',$CHR, $START, $END);
ok($slice->seq_region_name eq $CHR);
ok($slice->start == $START);
ok($slice->end   == $END);

#
# fetch_by_contig_name
#

my $projection = $slice->project('seqlevel');

#it is important to get a contig not cut off by slice start or end
unless(@$projection > 2) {
  warn("There aren't enough tiles in this path for this test to work");
}
my ($chr_start,$chr_end,$contig) = @{$projection->[1]};

ok($contig->length == ($chr_end - $chr_start + 1));

my $seq1 = $slice->subseq($chr_start, $chr_end);
my $seq2 = $contig->seq();

ok($seq1 eq $seq2);


#
# 12-13 fetch_by_fpc_name
#
#my $fpc_name = 'NT_011387';
#$slice = $slice_adaptor->fetch_by_supercontig_name($fpc_name);
#ok($new_slice->chr_start);
#ok($new_slice->chr_end);



#
# 14 - 15 fetch_by_clone_accession
#
#my $clone_acc = 'AL031658';
#$slice = $slice_adaptor->fetch_by_clone_accession($clone_acc);
#$new_slice = $slice_adaptor->fetch_by_clone_accession($clone_acc, $FLANKING);
#ok($new_slice->chr_start == $slice->chr_start - $FLANKING);
#ok($new_slice->chr_end   == $slice->chr_end   + $FLANKING);


#
# 16-17 fetch by transcript_stable_id
#
my $t_stable_id = 'ENST00000217315';
$slice = $slice_adaptor->fetch_by_transcript_stable_id($t_stable_id);
my $new_slice = $slice_adaptor->fetch_by_transcript_stable_id($t_stable_id,
                                                           $FLANKING);

ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);


#
# 18-19 fetch by transcript_id
#
my $transcript = $db->get_TranscriptAdaptor->fetch_by_stable_id($t_stable_id);
my $tid = $transcript->dbID;
$slice = $slice_adaptor->fetch_by_transcript_id($tid);
$new_slice = $slice_adaptor->fetch_by_transcript_id($tid, $FLANKING);
ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);


#
# 20-23 fetch_by_gene_stable_id
#
my $g_stable_id = 'ENSG00000125964';
$slice = $slice_adaptor->fetch_by_gene_stable_id($g_stable_id);
$new_slice = $slice_adaptor->fetch_by_gene_stable_id($g_stable_id, $FLANKING);
ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);

#verify we can retrieve the gene from this slice
my $gene_found = 0;
foreach my $g (@{$slice->get_all_Genes}) {
  if($g->stable_id eq $g->stable_id) {
    $gene_found = 1;
    last;
  }
}
ok($gene_found);

# same test for flanking slice
$gene_found = 0;
foreach my $g (@{$new_slice->get_all_Genes}) {
  if($g->stable_id eq $g->stable_id) {
    $gene_found = 1;
    last;
  }
}
ok($gene_found);


#
#  fetch_by_region (entire region)
#
$slice = $slice_adaptor->fetch_by_region('chromosome',$CHR);
ok($slice->seq_region_name eq $CHR);
ok($slice->start == 1);



#$slice = $slice_adaptor->fetch_by_chr_start_end("20", 29_252_000, 31_252_001 );
#my $name_list = $slice_adaptor->list_overlapping_supercontigs( $slice );

#for my $name ( @$name_list ) {
#  debug( "Overlapping supercontig ".$name );
#}

#ok( grep { $_ eq "NT_028392" } @$name_list);


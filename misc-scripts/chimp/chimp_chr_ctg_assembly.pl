use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(@ARGV);

my $slice_adaptor = $db->get_SliceAdaptor();

foreach my $slice (@{$slice_adaptor->fetch_all('chromosome')}) {
  my $asm_sr_id = $slice->get_seq_region_id();

  print STDERR ".";

  foreach my $segment (@{$slice->project('contig')}) {
    my $asm_start = $segment->from_start();
    my $asm_end   = $segment->from_end();
    my $cmp_sr_id = $segment->to_Slice->get_seq_region_id();
    my $cmp_start = $segment->to_Slice->start();
    my $cmp_end   = $segment->to_Slice->end();
    my $ori       = $segment->to_Slice->strand();

    print join("\t", $asm_sr_id, $cmp_sr_id, $asm_start, $asm_end,
               $cmp_start, $cmp_end, $ori), "\n";
  }
}

print STDERR "done\n";

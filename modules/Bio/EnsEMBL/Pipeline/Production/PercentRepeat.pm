package Bio::EnsEMBL::Pipeline::Production::PercentRepeat;

use base qw/Bio::EnsEMBL::Pipeline::Production::DensityGenerator/;
use Bio::EnsEMBL::Mapper::RangeRegistry;



use strict;
use warnings;


sub get_density {
  my ($self, $block, $slice) = @_;
  my $repeat_count = 0;
  my @repeats = @{ $block->get_all_RepeatFeatures() };
  @repeats = map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [$_->start, $_] } @repeats;
  my $curblock = undef;
  my @repeat_blocks;
  foreach my $repeat (@repeats) {
    if (defined($curblock) && $curblock->end >= $repeat->start) {
      if ($repeat->end > $block->length) {
        $curblock->end($block->length);
      } elsif ($repeat->end > $curblock->end) {
        $curblock->end($repeat->end);
      }
    } else {
      my $start = $repeat->start;
      my $end = $repeat->end;
      if ($repeat->end > $block->length) {
        $end = $block->length;
      }
      if ($repeat->start < 1) {
        $start = 1;
      }
      $curblock = Bio::EnsEMBL::Feature->new(-START => $start, -END => $end);
      push @repeat_blocks, $curblock;
    }
  }
  foreach my $repeat_block(@repeat_blocks) {
    $repeat_count += $repeat_block->length();
  }
  return 100*$repeat_count/$block->length();
}

return 1;


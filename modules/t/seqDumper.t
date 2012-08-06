use strict;
use warnings;

use Test::More;

use IO::String;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Utils::SeqDumper;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );


my $seq_dumper = Bio::EnsEMBL::Utils::SeqDumper->new();

ok(ref($seq_dumper) && $seq_dumper->isa('Bio::EnsEMBL::Utils::SeqDumper'));

my $file;

if($verbose) {
  $file = undef;
} else {
  $file = '/dev/null';
}


#do not dump snps they are not in core db
$seq_dumper->disable_feature_type('variation');

my $slice_adaptor = $db->get_SliceAdaptor();


my $slice =
  $slice_adaptor->fetch_by_region('contig', 'AL031658.11.1.162976');

$seq_dumper->dump($slice, 'EMBL', $file);
ok(1);

$seq_dumper->dump($slice, 'GENBANK', $file);
ok(1);

$seq_dumper->dump($slice, 'FASTA', $file);
ok(1);




$slice =
  $slice_adaptor->fetch_by_region('chromosome', '20', 30_500_000, 30_600_000);

$seq_dumper->dump($slice, 'EMBL', $file);
ok(1);

$seq_dumper->dump($slice, 'GENBANK', $file);
ok(1);

$seq_dumper->dump($slice, 'FASTA', $file);
ok(1);

my $index_fh = sub {
  my ($fh, $substr) = @_;
  $fh->setpos(0);
  my @lines;
  while(my $line = <$fh>) {
    chomp $line;
    push(@lines, $line) if index($line, $substr) == 0;
  }
  return \@lines;
};

my $index_count_fh = sub {
  my ($fh, $substr) = @_;
  return scalar(@{$index_fh->($fh, $substr)});
};

{
  my $frag_size = 1e7;
  my $seq = 'A'x$frag_size.'C'x$frag_size.'T'x$frag_size.'G'x$frag_size;
  my $sd = Bio::EnsEMBL::Utils::SeqDumper->new();
  $sd->{feature_types}->{$_} = 0 for keys %{$sd->{feature_types}};
  
  {
    my $fh = IO::String->new();
    $sd->dump_embl($slice, $fh, $seq);
    my $lines = $index_fh->($fh, 'SQ ');
    is(scalar(@{$lines}), 1, 'Expect only 1 EMBL SQ line describing a sequence');
    is($lines->[0], 'SQ   Sequence 40000000 BP; 10000000 A; 10000000 C; 10000000 G; 10000000 T; 0 other;', 'Formatting of SQ as expected');
  }
  
  {
    my $fh = IO::String->new();
    $sd->dump_genbank($slice, $fh, $seq);
    my $lines = $index_fh->($fh, 'BASE COUNT');
    is(@{$lines}, 1, 'Expect only 1 Genbank BASE COUNT line describing a sequence');
    is($lines->[0], 'BASE COUNT  10000000 a 10000000 c 10000000 g 10000000 t', 'Formatting of BASE COUNT as expected');
  }
  
}

done_testing();
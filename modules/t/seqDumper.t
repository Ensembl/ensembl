use strict;
use warnings;

use Bio::EnsEMBL::Utils::SeqDumper;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN { $| = 1;
	use Test;
	plan tests => 7;
}


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


use strict;
use warnings;

use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Utils::IO qw/:all/;

my $file = '/tmp/'.$ENV{USER}.'utilsIo.txt';
my $contents = <<'EOF';
>X
AAAAGGGTTCCC
TTGGCCAAAAAA
ATTC
EOF

throws_ok { slurp($file) } qr/No such file/, 'File does not currently exist so die';

work_with_file($file, 'w', sub {
  my ($fh) = @_;
  print $fh $contents;
  return;
});

my $written_contents = slurp($file);
is($contents, $written_contents, 'Contents should be the same');

my $written_contents_ref = slurp($file, 1);
is('SCALAR', ref($written_contents_ref), 'Asked for a ref so expect one back');
is($contents, $$written_contents_ref, 'Contents should be the same');

work_with_file($file, 'r', sub {
  my ($fh) = @_;
  my $line = <$fh>;
  chomp($line);
  is($line, '>X', 'First line expected to be FASTA header'); 
});

my $expected_array = [qw/>X AAAAGGGTTCCC TTGGCCAAAAAA ATTC/];
my $chomp = 1;
is_deeply(slurp_to_array($file, $chomp), $expected_array, 'Checking slurp to array with chomp');
$chomp = 0;
is_deeply(slurp_to_array($file, $chomp), [ map { "${_}\n" } @{$expected_array} ], 'Checking slurp to array with chomp');

unlink $file;

dies_ok { slurp($file) } 'File no longer exists so die';

done_testing();

use strict;
use warnings;

use Test::More;

BEGIN {
    use_ok('Bio::EnsEMBL::Utils::Iterator');
}

sub it_from_list {
    my @list = @_;
    return Bio::EnsEMBL::Utils::Iterator->new(sub{
        return shift @list;
    });
}

my ($it, $num, $i);

$it = it_from_list(1,2,3);

ok($it->has_next, 'iterator still has elements');

$i = 0;

while ($num = $it->next) {
    is($num, ++$i, 'got right element');
}

is($i, 3, 'got right number of elements');

# test empty iterator

my $empty_it = Bio::EnsEMBL::Utils::Iterator->new;

ok(!$empty_it->has_next, 'empty iterator has no elements');

# test grepping and mapping

$it = it_from_list(1,2,3,4);

my $even_it = $it->grep(sub {$_ % 2 == 0});

$i = 0;

while ($num = $even_it->next) {
    $i++;
    ok($num % 2 == 0, "$num is even");
}

is($i, 2, 'got right number of grepped elements back');

$it = it_from_list(1,2,3);

my $doubled_it = $it->map(sub {$_ * 2});

$i = 0;

while ($num = $doubled_it->next) {
    is($num, ++$i * 2, "$num = $i * 2");
}

is($i, 3, 'got right number of mapped elements back');

# test combining iterators

$it = it_from_list(1,2,3);

$empty_it = Bio::EnsEMBL::Utils::Iterator->new;

my $it2 = it_from_list(4,5,6);

my $combined_it = $it->append($empty_it, $it2);

$i = 0;

while ($num = $combined_it->next) {
    $i++;
}

is($i, 6, 'got right number of elements from combined iterators');

$it = it_from_list(1,2,3);

my $el = $it->peek;

is($el, $it->next, 'peek did not remove element');

# test converting iterator to an arrayref

$it = it_from_list(1,2,3);

is_deeply($it->to_arrayref, [1,2,3], 'got correct array back');

done_testing();


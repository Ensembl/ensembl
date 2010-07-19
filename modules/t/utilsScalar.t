use strict;
use warnings;

use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_ref);
use Bio::EnsEMBL::IdMapping::TinyGene;

my $gene = Bio::EnsEMBL::IdMapping::TinyGene->new_fast([]);

dies_ok { assert_ref(undef, 'ARRAY') } 'Undef value results in death';
dies_ok { assert_ref([], undef) } 'Undef assertion results in death';
throws_ok { assert_ref('string', 'ARRAY') } qr/produced no type/, 'Passing in a Scalar means death';
dies_ok { assert_ref(\'', 'ARRAY') } 'Ref of a Scalar is not an ARRAY so death';
dies_ok { assert_ref($gene, 'CODE') } 'TinyGene object is not a CODE so death';
dies_ok { assert_ref($gene, 'Bio::EnsEMBL::Feature') } 'TinyGene object is not a Bio::EnsEMBL::Feature so death';
dies_ok { assert_ref($gene, 'HASH') }  'TinyGene is blessed so we expect false even though it is a HASH';

lives_ok { assert_ref(\'', 'SCALAR') } 'Ref of a Scalar should be a SCALAR';
lives_ok { assert_ref([], 'ARRAY') } 'Ref of an array should be a ARRAY';
lives_ok { assert_ref({}, 'HASH') } 'Ref of a hash should be a HASH';
lives_ok { assert_ref($gene, 'Bio::EnsEMBL::IdMapping::TinyFeature') } 'Ref of a gene should be a TinyFeature';
lives_ok { assert_ref($gene, 'Bio::EnsEMBL::IdMapping::TinyGene') } 'Ref of a gene should be a TinyGene';

#Now for check_ref

dies_ok { check_ref([], undef) } 'Undef for assertion in check_ref results in death';

ok(! check_ref(undef, 'ARRAY'), 'Undef value returns false');
ok(! check_ref('string', 'ARRAY'), 'Passing in a Scalar means returns false');
ok(! check_ref(\'', 'ARRAY'),  'Ref of a Scalar is not an ARRAY so returns false');
ok(! check_ref($gene, 'CODE'),  'TinyGene object is not a CODE so returns false');
ok(! check_ref($gene, 'Bio::EnsEMBL::Feature'),  'TinyGene object is not a Bio::EnsEMBL::Feature so returns false');
ok(! check_ref($gene, 'HASH'),  'TinyGene is blessed so we expect false even though it is a HASH');

ok ( check_ref(\'', 'SCALAR'), 'Ref of a Scalar should be a SCALAR');
ok ( check_ref([], 'ARRAY'), 'Ref of an array should be a ARRAY');
ok ( check_ref({}, 'HASH'), 'Ref of a hash should be a HASH');
ok ( check_ref($gene, 'Bio::EnsEMBL::IdMapping::TinyFeature'), 'Ref of a gene should be a TinyFeature');
ok ( check_ref($gene, 'Bio::EnsEMBL::IdMapping::TinyGene'), 'Ref of a gene should be a TinyGene');

done_testing();

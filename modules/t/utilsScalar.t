# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Test::Exception;

use Bio::EnsEMBL::Utils::Scalar qw(:all);
use Bio::EnsEMBL::IdMapping::TinyGene;
use IO::Handle;

my $gene = Bio::EnsEMBL::IdMapping::TinyGene->new_fast([]);

# Assert ref check

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

# Array assertions
dies_ok { assert_array_contents([undef], 'ARRAY') } 'ARRAY: Undef value results in death';
dies_ok { assert_array_contents([], undef) } 'ARRAY: Undef assertion results in death';
throws_ok { assert_array_contents(['string'], 'ARRAY') } qr/produced no type/, 'ARRAY: Passing in a Scalar means death';
dies_ok { assert_array_contents([\''], 'ARRAY') } 'ARRAY: Ref of a Scalar is not an ARRAY so death';
dies_ok { assert_array_contents([$gene], 'CODE') } 'ARRAY: TinyGene object is not a CODE so death';
dies_ok { assert_array_contents([$gene], 'Bio::EnsEMBL::Feature') } 'ARRAY: TinyGene object is not a Bio::EnsEMBL::Feature so death';
dies_ok { assert_array_contents([$gene], 'HASH') }  'ARRAY: TinyGene is blessed so we expect false even though it is a HASH';

lives_ok { assert_array_contents([\''], 'SCALAR') } 'ARRAY: Ref of a Scalar should be a SCALAR';
lives_ok { assert_array_contents([[]], 'ARRAY') } 'ARRAY: Ref of an array should be a ARRAY';
lives_ok { assert_array_contents([{}], 'HASH') } 'ARRAY: Ref of a hash should be a HASH';
lives_ok { assert_array_contents([$gene], 'Bio::EnsEMBL::IdMapping::TinyFeature') } 'ARRAY: Ref of a gene should be a TinyFeature';
lives_ok { assert_array_contents([$gene], 'Bio::EnsEMBL::IdMapping::TinyGene') } 'ARRAY: Ref of a gene should be a TinyGene';
lives_ok { assert_array_contents([], 'Bio::EnsEMBL::IdMapping::TinyGene') } 'ARRAY: Empty array means no death';

# Array checks
dies_ok { check_array_contents([], undef) } 'ARRAY: Undef for assertion in check_array_contents results in death';

ok(! check_array_contents([undef], 'ARRAY'), 'ARRAY: Undef value returns false');
ok(! check_array_contents(['string'], 'ARRAY'), 'ARRAY: Passing in a Scalar means returns false');
ok(! check_array_contents([\''], 'ARRAY'),  'ARRAY: Ref of a Scalar is not an ARRAY so returns false');
ok(! check_array_contents([$gene], 'CODE'),  'ARRAY: TinyGene object is not a CODE so returns false');
ok(! check_array_contents([$gene], 'Bio::EnsEMBL::Feature'),  'ARRAY: TinyGene object is not a Bio::EnsEMBL::Feature so returns false');
ok(! check_array_contents([$gene], 'HASH'),  'ARRAY: TinyGene is blessed so we expect false even though it is a HASH');

ok ( check_array_contents([\''], 'SCALAR'), 'ARRAY: Ref of a Scalar should be a SCALAR');
ok ( check_array_contents([[]], 'ARRAY'), 'ARRAY: Ref of an array should be a ARRAY');
ok ( check_array_contents([{}], 'HASH'), 'ARRAY: Ref of a hash should be a HASH');
ok ( check_array_contents([$gene], 'Bio::EnsEMBL::IdMapping::TinyFeature'), 'ARRAY: Ref of a gene should be a TinyFeature');
ok ( check_array_contents([$gene], 'Bio::EnsEMBL::IdMapping::TinyGene'), 'ARRAY: Ref of a gene should be a TinyGene');

# Hash assertions
dies_ok { assert_hash_contents({a => undef}, 'ARRAY') } 'HASH: Undef value results in death';
dies_ok { assert_hash_contents({}, undef) } 'HASH: Undef assertion results in death';
throws_ok { assert_hash_contents({ a => 'string'}, 'ARRAY') } qr/produced no type/, 'HASH: Passing in a Scalar means death';
dies_ok { assert_hash_contents({a => \''}, 'ARRAY') } 'HASH: Ref of a Scalar is not an ARRAY so death';
dies_ok { assert_hash_contents({a => $gene}, 'CODE') } 'HASH: TinyGene object is not a CODE so death';
dies_ok { assert_hash_contents({a => $gene}, 'Bio::EnsEMBL::Feature') } 'HASH: TinyGene object is not a Bio::EnsEMBL::Feature so death';
dies_ok { assert_hash_contents({a => $gene}, 'HASH') }  'HASH: TinyGene is blessed so we expect false even though it is a HASH';

lives_ok { assert_hash_contents({a => \'' }, 'SCALAR') } 'HASH: Ref of a Scalar should be a SCALAR';
lives_ok { assert_hash_contents({a => []}, 'ARRAY') } 'HASH: Ref of an array should be a ARRAY';
lives_ok { assert_hash_contents({a => {}}, 'HASH') } 'HASH: Ref of a hash should be a HASH';
lives_ok { assert_hash_contents({a => $gene}, 'Bio::EnsEMBL::IdMapping::TinyFeature') } 'HASH: Ref of a gene should be a TinyFeature';
lives_ok { assert_hash_contents({a => $gene}, 'Bio::EnsEMBL::IdMapping::TinyGene') } 'HASH: Ref of a gene should be a TinyGene';
lives_ok { assert_hash_contents({}, 'Bio::EnsEMBL::IdMapping::TinyGene') } 'HASH: Empty array means no death';

# Hash checks
dies_ok { check_hash_contents({}, undef) } 'HASH: Undef for assertion in check_hash_contents results in death';

ok(! check_hash_contents({a => undef}, 'ARRAY'), 'HASH: Undef value returns false');
ok(! check_hash_contents({a => 'string'}, 'ARRAY'), 'HASH: Passing in a Scalar means returns false');
ok(! check_hash_contents({a => \''}, 'ARRAY'),  'HASH: Ref of a Scalar is not an ARRAY so returns false');
ok(! check_hash_contents({a => $gene}, 'CODE'),  'HASH: TinyGene object is not a CODE so returns false');
ok(! check_hash_contents({a => $gene}, 'Bio::EnsEMBL::Feature'),  'HASH: TinyGene object is not a Bio::EnsEMBL::Feature so returns false');
ok(! check_hash_contents({a => $gene}, 'HASH'),  'HASH: TinyGene is blessed so we expect false even though it is a HASH');

ok ( check_hash_contents({a => \''}, 'SCALAR'), 'HASH: Ref of a Scalar should be a SCALAR');
ok ( check_hash_contents({a => []}, 'ARRAY'), 'HASH: Ref of an array should be a ARRAY');
ok ( check_hash_contents({a => {}}, 'HASH'), 'HASH: Ref of a hash should be a HASH');
ok ( check_hash_contents({a => $gene}, 'Bio::EnsEMBL::IdMapping::TinyFeature'), 'HASH: Ref of a gene should be a TinyFeature');
ok ( check_hash_contents({a => $gene}, 'Bio::EnsEMBL::IdMapping::TinyGene'), 'HASH: Ref of a gene should be a TinyGene');


#Array Wrapping

my $undef_ref = undef;
my $value = 'hello';
my $array = [$value];
is_deeply( [], wrap_array(), 'Checking empty value means empty array');
is_deeply( [], wrap_array(undef), 'Checking undef means empty array');
is_deeply( [], wrap_array($undef_ref), 'Checking undef ref means empty array');
is_deeply( [$value], wrap_array($value), 'Checking Scalar ref means wrapped array');
is_deeply( $array, wrap_array($array), 'Checking arrays are the same if given array ref');
is( $array, wrap_array($array), 'Checking arrays are the same reference');
is_deeply( [{a => $value}], wrap_array({ a => $value}), 'Checking code behaves when working with hashes');

#Ref Can
my $blessed_array = bless([], 'Bio::EnsEMBL::BrianBlessedArray');
throws_ok { check_ref_can('string', undef) } qr/method/, 'Passing in no method means death';
ok(! check_ref_can(undef, 'met'), 'Passing in an undefined value means false');
ok(! check_ref_can('string', 'met'), 'Passing in an unblessed value means false');
ok(check_ref_can($gene, 'start'), 'TinyGene implements start()');
ok(!check_ref_can($gene, 'wibble'), 'TinyGene does not implement wibble()');
ok(!check_ref_can($blessed_array, 'wibble'), 'The blessed array does not implement any methods let alone wibble()');

#Ref assert can
throws_ok { assert_ref_can('string', undef) } qr/method/, 'Passing in no method means death';
dies_ok { assert_ref_can(undef, 'met')} 'Passing in an undefined value means death';
dies_ok { assert_ref_can('string', 'met')} 'Passing in an unblessed value means death';
lives_ok { assert_ref_can($gene, 'start')} 'TinyGene implements start()';
dies_ok { assert_ref_can($gene, 'wibble')} 'TinyGene does not implement wibble() so death';
dies_ok { assert_ref_can($blessed_array, 'wibble')} 'The blessed array does not implement any methods let alone wibble() so death';

#Numerics
throws_ok { assert_numeric(undef) } qr/undefined/, 'Passing in undefined scalar means death';
dies_ok { assert_numeric(bless(1, 'Brian'), 'met')} 'Passing in a blessed scalar means death';
dies_ok { assert_numeric('hello')} 'Passing in a String scalar means death';
dies_ok { assert_numeric({})} 'Passing in a HashRef means death';
lives_ok { assert_numeric(1E-10) } 'Passing in scientific notation numeric means lives';
lives_ok { assert_numeric(1.2) } 'Passing in floating point means lives';
lives_ok { assert_numeric(1) } 'Passing in integer means lives';

#Integers
throws_ok { assert_integer(undef) } qr/undefined/, 'Passing in undefined scalar means death';
dies_ok { assert_integer(bless(1, 'Brian'), 'met')} 'Passing in a blessed scalar means death';
dies_ok { assert_integer('hello')} 'Passing in a String scalar means death';
dies_ok { assert_integer({})} 'Passing in a HashRef means death';
dies_ok { assert_integer(1E-10) } 'Passing in negative scientific notation numeric means death';
lives_ok { assert_integer(1E10) } 'Passing in positive scientific notation numeric means lives';
lives_ok { assert_integer(1_000) } 'Separators means lives';
lives_ok { assert_integer('2') } 'A string numeric is ok';
dies_ok { assert_integer(1.2) } 'Passing in floating point means death';
lives_ok { assert_integer(1) } 'Passing in integer means lives';

#Strand
throws_ok { assert_strand(undef) } qr/undefined/, 'Passing in undefined scalar means death';
dies_ok { assert_strand(bless(1, 'Brian'), 'met')} 'Passing in a blessed scalar means death';
dies_ok { assert_strand('hello')} 'Passing in a String scalar means death';
dies_ok { assert_strand({})} 'Passing in a HashRef means death';
dies_ok { assert_strand(1E-10) } 'Passing in scientific notation numeric means death';
dies_ok { assert_strand(1.2) } 'Passing in floating point means death';
dies_ok { assert_strand(2) } 'Passing in floating point means death';
lives_ok { assert_strand(1) } 'Passing in integer 1 means lives';
lives_ok { assert_strand(0) } 'Passing in integer 0 means lives';
lives_ok { assert_strand(-1) } 'Passing in integer -1 means lives';

#Boolean
throws_ok { assert_boolean(undef) } qr/undefined/, 'Passing in undefined scalar means death';
dies_ok { assert_boolean(bless(1, 'Brian'), 'met')} 'Passing in a blessed scalar means death';
dies_ok { assert_boolean('hello')} 'Passing in a String scalar means death';
dies_ok { assert_boolean({})} 'Passing in a HashRef means death';
dies_ok { assert_boolean(1E-10) } 'Passing in scientific notation numeric means death';
dies_ok { assert_boolean(1.2) } 'Passing in floating point means death';
dies_ok { assert_boolean(-1) } 'Passing in integer -1 means death';
lives_ok { assert_strand(1) } 'Passing in integer 1 means lives';
lives_ok { assert_strand(0) } 'Passing in integer 0 means lives';

#File handles
my $scalar;
my $other_scalar;
open my $scalar_fh, '>', \$scalar;
open my $other_scalar_fh, '>', \$other_scalar;
bless($other_scalar_fh, __PACKAGE__);
my $io_handle = IO::Handle->new(); # no need to close as it isn't opened yet just created
throws_ok { assert_file_handle(undef) } qr/undefined/, 'Passing in undefined scalar means death';
dies_ok { assert_file_handle(bless(1, 'Brian'), 'met')} 'Passing in a blessed scalar means death';
dies_ok { assert_file_handle('hello')} 'Passing in a String scalar means death';
dies_ok { assert_file_handle({})} 'Passing in a HashRef means death';
dies_ok { assert_file_handle(1E-10) } 'Passing in scientific notation numeric means death';
dies_ok { assert_file_handle(1.2) } 'Passing in floating point means death';
dies_ok { assert_file_handle(-1) } 'Passing in integer -1 means death';
dies_ok { assert_file_handle(1) } 'Passing in integer 1 means death';
lives_ok { assert_file_handle($scalar_fh) } 'Passing in a scalar FH means lives';
lives_ok { assert_file_handle($other_scalar_fh) } 'Passing in a blessed scalar FH means lives';
lives_ok { assert_file_handle($io_handle) } 'Passing in an IO::Handle means lives';
close($_) for ($scalar_fh, $other_scalar_fh);

#Scope Guard

#First normal circumstances
{
  my $v = 'wibble';
  is($v, 'wibble', 'Value is normal');
  {
    my $guard = scope_guard(sub { $v = 'wibble'});
    $v = 'wobble';
    is($v, 'wobble', 'Value has been changed');
  }
  is($v, 'wibble', 'Value has been reset');
}

#With a die; Perl respects DESTROY even if we are in the middle of SIG{__DIE__}
{
  my $v = 'wibble';
  is($v, 'wibble', 'Value is normal');
  $v = 'wobble';
  is($v, 'wobble', 'Value has been changed');
  eval {
    my $guard = scope_guard(sub { $v = 'wibble'});
    die 'we have died';
  };
  is($v, 'wibble', 'Value has been reset even after a die');
}

#Array split
{
  my $original_array = [1..7];
  my $split_two = split_array(2, $original_array);
  is_deeply($split_two, [[1,2],[3,4],[5,6],[7]], 'Checking split of 7 element array into arrays of max size 2') or diag explain $split_two;
  my $split_ten = split_array(10, $original_array);
  is_deeply($split_ten, [$original_array], 'Checking split of 7 element array into arrays of max size 10') or diag explain $split_ten;
  dies_ok { split_array(1, {}) } 'Passing in a non-array ref means death';
}

# Testing global assertion turn-off
{
  local $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0;
  lives_ok {assert_ref([], 'HASH')      } 'Assertions off; [] returns true for HASH';
  lives_ok {assert_file_handle([])      } 'Assertions off; [] returns true for assert_file_handle()';
  lives_ok {assert_strand([])           } 'Assertions off; [] returns true for assert_strand()';
  lives_ok {assert_boolean([])          } 'Assertions off; [] returns true for assert_boolean()';
  lives_ok {assert_integer([])          } 'Assertions off; [] returns true for assert_integer()';
  lives_ok {assert_numeric([])          } 'Assertions off; [] returns true for assert_numeric()';
  lives_ok {assert_ref_can([], 'wibble')} 'Assertions off; [] returns true for assert_ref_can()';
  lives_ok {assert_array_contents([{}], 'wibble')} 'Assertions off; [{}] returns true for assert_array_contents() when asserting type "wibble"';
  lives_ok {assert_hash_contents({a => {}}, 'wibble')} 'Assertions off; {a => {}} returns true for assert_hash_contents() when asserting type "wibble"';
}
dies_ok { assert_ref([], 'HASH') } 'Assertions back on; [] is not a HASH';

done_testing();

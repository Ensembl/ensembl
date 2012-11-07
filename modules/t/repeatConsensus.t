use strict;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::RepeatConsensus;

use Test::More;
my $verbose = 0;

my $consensus = 'actg';
my $name   =  'actg(n)';
my $length =  4;
my $class  = 'Simple_repeat';
my $dbID  = 123;

#
# Test constructor
#
my $repeat_consensus = Bio::EnsEMBL::RepeatConsensus->new
  (-REPEAT_CONSENSUS => $consensus,
   -NAME             => $name,
   -LENGTH           => $length,
   -REPEAT_CLASS    => $class,
   -DBID             => 123);

ok ($repeat_consensus && ref($repeat_consensus) && 
    $repeat_consensus->isa('Bio::EnsEMBL::RepeatConsensus'));

ok($repeat_consensus->length() == $length);
ok($repeat_consensus->repeat_consensus() eq $consensus);
ok($repeat_consensus->seq() eq $consensus);
ok($repeat_consensus->name() eq $name);
ok($repeat_consensus->dbID() == $dbID);
ok($repeat_consensus->repeat_class() eq $class);

ok(test_getter_setter($repeat_consensus,'length',10));
ok(test_getter_setter($repeat_consensus,'repeat_class','dummy'));
ok(test_getter_setter($repeat_consensus,'name','dummy'));
ok(test_getter_setter($repeat_consensus,'repeat_consensus','ATGCATGCAT'));
ok(test_getter_setter($repeat_consensus,'dbID',42));

done_testing();
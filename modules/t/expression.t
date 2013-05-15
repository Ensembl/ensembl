
use strict;
use warnings;

use Bio::EnsEMBL::Expression;

use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts

use Test::More;

use Bio::EnsEMBL::Test::TestUtils;

#
# test constructor
#

my $name     = 'testname';
my $desc     = 'testdesc';
my $ontology = 'testontology';
my $value    = 'testval';

my $tissue = Bio::EnsEMBL::Expression->new
  (-NAME => $name,
   -DESCRIPTION => $desc,
   -ONTOLOGY => $ontology,
   -VALUE => $value);

ok($tissue->name()  eq $name);
ok($tissue->description()  eq $desc);
ok($tissue->ontology()  eq $ontology);
ok($tissue->value() eq $value);

#
# test getter/setters
#
ok(test_getter_setter($tissue, 'name', 'newname'));
ok(test_getter_setter($tissue, 'ontology', 'newontology'));
ok(test_getter_setter($tissue, 'description', 'newdesc'));
ok(test_getter_setter($tissue, 'value', 'newvalue'));

done_testing();


use strict;
use warnings;

use Bio::EnsEMBL::Attribute;

use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts

BEGIN { $| = 1;
	use Test;
	plan tests => 8;
}

use Bio::EnsEMBL::Test::TestUtils;

#
# test constructor
#

my $code = 'testcode';
my $name = 'testname';
my $desc = 'testdesc';
my $value = 'testval';

my $attrib = Bio::EnsEMBL::Attribute->new
  (-CODE => $code,
   -NAME => $name,
   -DESCRIPTION => $desc,
   -VALUE => $value);

ok($attrib->code()  eq $code);
ok($attrib->name()  eq $name);
ok($attrib->description()  eq $desc);
ok($attrib->value() eq $value);

#
# test getter/setters
#
ok(test_getter_setter($attrib, 'name', 'newname'));
ok(test_getter_setter($attrib, 'code', 'newcode'));
ok(test_getter_setter($attrib, 'description', 'newdesc'));
ok(test_getter_setter($attrib, 'value', 'newvalue'));


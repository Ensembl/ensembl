use lib 't';

use strict;
use warnings;

use Bio::EnsEMBL::Attribute;

use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts

BEGIN { $| = 1;
	use Test;
	plan tests => 8;
}

use TestUtils qw( debug test_getter_setter );

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
   -DESC => $desc,
   -VALUE => $value);

ok($attrib->code()  eq $code);
ok($attrib->name()  eq $name);
ok($attrib->desc()  eq $desc);
ok($attrib->value() eq $value);

#
# test getter/setters
#
ok(test_getter_setter($attrib, 'name', 'newname'));
ok(test_getter_setter($attrib, 'code', 'newcode'));
ok(test_getter_setter($attrib, 'desc', 'newdesc'));
ok(test_getter_setter($attrib, 'value', 'newvalue'));


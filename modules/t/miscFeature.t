use lib 't';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 11
}

use MultiTestDB;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::MiscSet;
use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts


#test constructor
my $mf = Bio::EnsEMBL::MiscFeature->new(-START => 10,
                                        -END   => 100);

ok($mf->start() == 10 && $mf->end() == 100);



#
# Test add_set, get_set, get_set_codes
#
my $ms1 = Bio::EnsEMBL::MiscSet->new(3, undef,
                                     '1mbcloneset',
                                     '1mb Clone Set',
                                     'This is a 1MB cloneset',
                                     1e7);

my $ms2 = Bio::EnsEMBL::MiscSet->new(4, undef,
                                     'tilepath',
                                     'Tiling Path',
                                     'NCBI33 Tiling Path',
                                     1e6);



$mf->add_set($ms1);
$mf->add_set($ms2);


ok($mf->get_set($ms1->code) == $ms1);
ok($mf->get_set($ms2->code) == $ms2);

my @codes = $mf->get_set_codes;

ok(@codes == 2);

@codes = grep {($_ eq $ms1->code()) || ($_ eq $ms2->code())} @codes;

ok(@codes == 2);

#
# Test add_attribute, get_attribute_types, get_attribute
#

my $name1 = 'test name';
my $name2 = 'AL4231124.1';

$mf->add_attribute('name',  $name1);

ok($mf->display_id eq $name1);

$mf->add_attribute('name',  $name2);
$mf->add_attribute('version', '4');

my @types = $mf->get_attribute_types();
ok(@types == 2);

@types = grep {$_ eq 'name' || $_ eq 'version'} @types;

ok(@types == 2);

my @attribs = $mf->get_attribute('name');

ok(@attribs == 2);

@attribs = grep {$_ eq $name1 || $_ eq $name2} @attribs;

ok(@attribs == 2);

@attribs = $mf->get_attribute('version');
ok(@attribs == 1 && $attribs[0] eq '4');

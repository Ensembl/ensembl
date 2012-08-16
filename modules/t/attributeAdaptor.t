use strict;

use Bio::EnsEMBL::Test::TestUtils;

use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Attribute;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $db = $multi->get_DBAdaptor("core");

my $slice_adaptor = $db->get_SliceAdaptor();
my $mfa           = $db->get_MiscFeatureAdaptor();


#
# Test get_AttributeAdaptor works
#
my $aa = $db->get_AttributeAdaptor();

ok($aa && ref($aa) && $aa->isa('Bio::EnsEMBL::DBSQL::AttributeAdaptor'));


# hide the contents of the attrib_type, misc_attrib, seq_region_attrib tables
# so we can test storing etc. with a clean slate
$multi->hide('core', 'misc_attrib', 'seq_region_attrib', 'attrib_type');


##############
# MiscFeature functionality tests
#

my $attrib = Bio::EnsEMBL::Attribute->new(-NAME => 'test_name',
                                          -CODE => 'test_code',
                                          -DESCRIPTION => 'test_desc',
                                          -VALUE => 'test_value');


my $mf = $mfa->fetch_by_dbID(1);
$aa->store_on_MiscFeature($mf, [$attrib]);

#
# make sure the misc_attrib table was updated
#
my $count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM misc_attrib " .
   "WHERE misc_feature_id = 1")->[0]->[0];

ok($count == 1);

#
# make sure the attrib_type table was updated
#
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM attrib_type " .
   "WHERE code = 'test_code'")->[0]->[0];
ok($count == 1);

#
# test that we can now retrieve this attribute
#
my @attribs = @{$aa->fetch_all_by_MiscFeature($mf)};

ok(@attribs == 1);

$attrib = $attribs[0];

ok($attrib->name eq 'test_name');
ok($attrib->code eq 'test_code');
ok($attrib->description eq 'test_desc');
ok($attrib->value eq 'test_value');


#
# test the removal of this attribute
#
$aa->remove_from_MiscFeature($mf);

#
# make sure the misc_attrib table was updated
#
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM misc_attrib " .
   "WHERE misc_feature_id = 1")->[0]->[0];

ok($count == 0);

#
# make sure the attribute is no longer retrievable
#
@attribs = @{$aa->fetch_all_by_MiscFeature($mf)};
ok(@attribs == 0);



#################
# Slice functionality tests
#


$attrib = Bio::EnsEMBL::Attribute->new(-NAME => 'test_name2',
                                       -CODE => 'test_code2',
                                       -DESCRIPTION => 'test_desc2',
                                       -VALUE => 'test_value2');


my $slice = $slice_adaptor->fetch_by_region('chromosome', '20');

$aa->store_on_Slice($slice, [$attrib]);

#
# make sure the seq_region_attrib table was updated
#
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM seq_region_attrib " .
   "WHERE seq_region_id = ".$slice->get_seq_region_id())->[0]->[0];

ok($count == 1);

#
# make sure the attrib_type table was updated
#
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM attrib_type " .
   "WHERE code = 'test_code2'")->[0]->[0];
ok($count == 1);

#
# test that we can now retrieve this attribute
#
@attribs = @{$aa->fetch_all_by_Slice($slice)};

ok(@attribs == 1);


@attribs = @{$aa->fetch_all_by_Slice($slice,"rubbish")};
ok(@attribs == 0);

@attribs = @{$aa->fetch_all_by_Slice($slice,"test_code2")};
ok(@attribs == 1);



$attrib = $attribs[0];

ok($attrib->name eq 'test_name2');
ok($attrib->code eq 'test_code2');
ok($attrib->description eq 'test_desc2');
ok($attrib->value eq 'test_value2');


#
# test the removal of this attribute with atrrib code
#
$aa->remove_from_Slice($slice,"junk");
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM seq_region_attrib " .
   "WHERE seq_region_id = " . $slice->get_seq_region_id())->[0]->[0];

ok($count == 1);




#
# test the removal of this attribute
#

$aa->remove_from_Slice($slice,"test_code2");
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM seq_region_attrib " .
   "WHERE seq_region_id = " . $slice->get_seq_region_id())->[0]->[0];

ok($count == 0);



#
# make sure the attribute is no longer retrievable
#
@attribs = @{$aa->fetch_all_by_Slice($slice)};
ok(@attribs == 0);



#
# try to add an attribute with an already existing code
#
$aa->store_on_Slice($slice, [$attrib]);
#
# make sure the seq_region_attrib table was updated
#
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM seq_region_attrib " .
   "WHERE seq_region_id = ".$slice->get_seq_region_id())->[0]->[0];

ok($count == 1);

#
# make sure the attrib_type table was updated
#
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM attrib_type " .
   "WHERE code = 'test_code2'")->[0]->[0];
ok($count == 1);

@attribs = @{$aa->fetch_all_by_Slice($slice)};
note "attribs: " . scalar(@attribs);
ok(@attribs == 1);


#
# test the removal of this attribute
#
$aa->remove_from_Slice($slice);
$count = $db->dbc->db_handle->selectall_arrayref
  ("SELECT count(*) FROM seq_region_attrib " .
   "WHERE seq_region_id = " . $slice->get_seq_region_id())->[0]->[0];

ok($count == 0);

#
# test the storage of empty attrib values
#
{
  my %args = (-NAME => 'test_name2', -CODE => 'test_code2', -DESCRIPTION => 'test_desc2');
  my $current_rows = count_rows($db, 'seq_region_attrib');
  my $atrib = Bio::EnsEMBL::Attribute->new(%args,);
  $aa->store_on_Slice($slice, [Bio::EnsEMBL::Attribute->new(%args, -VALUE => q{})]);
  $aa->store_on_Slice($slice, [Bio::EnsEMBL::Attribute->new(%args, -VALUE => 0)]);
  my $new_rows = count_rows($db, 'seq_region_attrib');
  cmp_ok($new_rows, '>', $current_rows, 'Asserting the storage of undefined attributes will always store them');
}

$multi->restore('core', 'misc_attrib', 'seq_region_attrib', 'attrib_type');

done_testing();

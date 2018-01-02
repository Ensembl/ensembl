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

use Bio::EnsEMBL::Test::TestUtils;

use Test::More;
use Test::Warnings;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::MiscSet;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $dba = $multi->get_DBAdaptor("core");


#
# Test get_MiscFeatureAdaptor works
#
my $mfa = $dba->get_MiscFeatureAdaptor();

ok($mfa && ref($mfa) && $mfa->isa('Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor'));


#
# Test fetching by slice
#


my $chr_slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome', '20');
my $features = $mfa->fetch_all_by_Slice($chr_slice);
debug('--- chr 20 misc_features ---');
debug("Got " . scalar(@$features));
ok(@$features == 7);
print_features($features);


$features = $mfa->fetch_all_by_Slice_and_set_code($chr_slice,'ntctgs');
debug('--- chr 20 ntcontigs set---');
debug("Got " . scalar(@$features));
ok(@$features == 3);
print_features($features);

$features = $mfa->fetch_all_by_Slice_and_set_code($chr_slice,'cloneset');
debug('--- chr 20 cloneset set---');
debug("Got " . scalar(@$features));
ok(@$features == 0);
print_features($features);

my $feature = $mfa->fetch_by_dbID(1);
debug('--- fetch_by_dbID ---');
ok($feature->dbID() == 1);
ok($feature->start() == 61140848);
ok($feature->end()   == 62842997);
ok($feature->strand() == 1);
print_features([$feature]);


#
# Test fetching by attribute
#

debug("--- fetch by attribute (superctg) ---");
$features = $mfa->fetch_all_by_attribute_type_value('superctg');

ok(@$features == 7);
print_features($features);

debug("--- fetch by attribute (superctg, NT_035608) ---");
$features = $mfa->fetch_all_by_attribute_type_value('superctg','NT_035608');

ok(@$features == 1);
print_features($features);

debug("--- fetch by attribute (embl_acc) ---");
$features = $mfa->fetch_all_by_attribute_type_value('embl_acc');
ok(@$features == 1);
print_features($features);

#
# Test fetch_by_attribute_set_value method
#
my $feature = $mfa->fetch_by_attribute_set_value('embl_acc', 'AL000001', 'ntctgs');
debug('--- fetch_by_attribute_set_value ---');
ok($feature->dbID() == 2);
print_features([$feature]);

$feature = $mfa->fetch_by_attribute_set_value('embl_acc', 'rubbish', 'ntctgs');
ok(!defined $feature);

$feature = $mfa->fetch_by_attribute_set_value('embl_acc', 'AL000001', 'rubbish');
ok(!defined $feature);

$feature = $mfa->fetch_by_attribute_set_value('rubbish', 'AL000001', 'ntctgs');
ok(!defined $feature);


$multi->hide('core', 'misc_feature', 'misc_feature_misc_set', 'meta_coord',
            'misc_attrib', 'attrib_type', 'misc_set');

#
# test store method
#
my $misc_set = Bio::EnsEMBL::MiscSet->new
  (-NAME => 'setname',
   -CODE => 'setcode',
   -DESCRIPTION => 'setdescription',
   -LONGEST_FEATURE => 10000);

my $attrib = Bio::EnsEMBL::Attribute->new
  (-NAME => 'attribute',
   -CODE => 'attribcode',
   -DESCRIPTION => 'attribdescription',
   -VALUE => 'testvalue');

my $mf = Bio::EnsEMBL::MiscFeature->new
  (-START  => 100,
   -END    => 200,
   -STRAND => 1,
   -SLICE  => $chr_slice);


$mf->add_MiscSet($misc_set);
$mf->add_Attribute($attrib);
$mfa->store($mf);

ok($mf->dbID() && $mf->adaptor());
ok($misc_set->dbID() && $misc_set->adaptor());

$mf = $mfa->fetch_by_dbID($mf->dbID());

my @attribs = @{$mf->get_all_Attributes()};
ok(@attribs == 1);
ok($attribs[0]->code eq 'attribcode');

my @sets = @{$mf->get_all_MiscSets()};
ok(@sets == 1);
ok($sets[0]->code eq 'setcode');



# try to store a misc feature without attributes

$mf = Bio::EnsEMBL::MiscFeature->new
  (-START => 100,
   -END  => 200,
   -STRAND => 1,
   -SLICE => $chr_slice);

$mfa->store($mf);

ok($mf->is_stored($dba));


$multi->restore('core', 'misc_feature', 'misc_feature_misc_set', 'meta_coord',
                'misc_attrib', 'attrib_type', 'misc_set');





sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      my @attribs = @{$f->get_all_Attributes() };
      my @sets = @{$f->get_all_MiscSets()};

      
      my $attrib_string = join(":", map { $_->code()." => ".$_->value() } @attribs );;
      my $set_string = join(":", map { $_->code() } @sets );

      my $seqname = $f->slice->seq_region_name();
      debug($seqname . ' ' . $f->start().'-'.$f->end().'('.$f->strand().
            ') ['. $f->dbID.'] '." ($attrib_string) ($set_string)");
    } else {
      debug('UNDEF');
    }
  }
}

done_testing();

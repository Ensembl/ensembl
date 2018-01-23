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
use Bio::EnsEMBL::SimpleFeature;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
 
# get a core DBAdaptor
#
my $dba = $multi->get_DBAdaptor("core");
my $sfa = $dba->get_SimpleFeatureAdaptor;

#
# 1 create a new Simplefeature
#
my $sf = new Bio::EnsEMBL::SimpleFeature;
ok($sf);


#
# 2-7 test the basic getter and setters
#

# 2 start
ok(test_getter_setter($sf,'start',10));

# 3 end
ok(test_getter_setter($sf,'end',14));

# 4 strand
ok(test_getter_setter($sf,'strand',1));

# 5 score
ok(test_getter_setter($sf,'score',42));

# 6 display_label
ok(test_getter_setter($sf,'display_label','dummy_label'));

# 7 dbID
ok(test_getter_setter($sf,'dbID',42));



#
# 9 check adaptor attaching
#
$sf->adaptor($sfa);
ok($sf->adaptor->isa('Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor'));


my $chr_slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome', '20');
my $features = $sfa->fetch_all_by_Slice($chr_slice);
debug('--- chr 20 simple features ---');
debug("Got " . scalar(@$features));
ok(@$features == 136);
print_features($features);


my $cln_slice = $dba->get_SliceAdaptor->fetch_by_region('clone','AL031658.11');
$features = $sfa->fetch_all_by_Slice($cln_slice);
debug('-- cln AL031658.11 simple features ---');
debug("Got " . scalar(@$features));
ok(@$features == 20);
print_features($features);

my $sprctg_slice = $dba->get_SliceAdaptor->fetch_by_region('supercontig',
                                                         'NT_028392');
$features = $sfa->fetch_all_by_Slice($sprctg_slice);
debug('-- sprctg NT_028392 simple features ---');
debug("Got " . scalar(@$features));
ok(@$features == 136);
print_features($features);

my $ctg_slice = $dba->get_SliceAdaptor->fetch_by_region('contig',
                                                       'AL031658.11.1.162976');
$features = $sfa->fetch_all_by_Slice($ctg_slice);
debug('--- contig AL031658.11.1.162976 simple features ---');
debug("Got " . scalar(@$features));
ok(@$features == 20);
print_features($features);



#retrieve a feature via dbID

debug('---- fetch_by_dbID (default coords) ----');
my $feat = $sfa->fetch_by_dbID(14564);
ok($feat->dbID == 14564);
ok($feat->slice->seq_region_name() eq 'AL031658.11.1.162976');
ok($feat->start == 64109);
ok($feat->end == 64112);
ok($feat->strand == 1);

print_features([$feat]);

#transform to chromosome coords
debug('---- tranform to chromosomal ----');
$feat = $feat->transform('chromosome');
ok($feat->dbID == 14564);
ok($feat->slice->seq_region_name() eq '20');
ok($feat->start == 30327623);
ok($feat->end == 30327626);
ok($feat->strand == 1);
print_features([$feat]);

# transform to supercontig
debug('---- tranform to supercontig ----');
$feat = $feat->transform('contig')->transform('supercontig');
ok($feat->dbID == 14564);
ok($feat->slice->seq_region_name() eq 'NT_028392');
ok($feat->start == 735658);
ok($feat->end == 735661);
ok($feat->strand == 1);
print_features([$feat]);

debug('---- transform to clone ----');
$feat = $feat->transform('contig')->transform('clone');
ok($feat->dbID == 14564);
ok($feat->slice->seq_region_name() eq 'AL031658.11');
ok($feat->start == 64109);
ok($feat->end == 64112);
ok($feat->strand == 1);
print_features([$feat]);


# List_dbidx
my $ids = $sfa->list_dbIDs();
ok (@{$ids});


ok($feat->display_id eq $feat->display_label);



sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      my $seqname = $f->slice->seq_region_name();
      my $analysis = $f->analysis->logic_name();
      debug($seqname . ' ' . $f->start().'-'.$f->end().'('.$f->strand().
            ') ['. $f->dbID.'] '. $f->display_label.' '.$f->score() .
            " ($analysis)");
    } else {
      debug('UNDEF');
    }
  }
}

done_testing();

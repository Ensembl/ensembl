# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::SplicingEvent;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
 
# get a core DBAdaptor
#
my $dba = $multi->get_DBAdaptor("patch");
my $sfa = $dba->get_SplicingEventAdaptor;

#
# 1 create a new Splicingevent
#
my $sf = new Bio::EnsEMBL::SplicingEvent;
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

# 6 display_label
ok(test_getter_setter($sf,'name','dummy_label'));

# 7 dbID
ok(test_getter_setter($sf,'dbID',42));



#
# 9 check adaptor attaching
#
$sf->adaptor($sfa);
ok($sf->adaptor->isa('Bio::EnsEMBL::DBSQL::SplicingEventAdaptor'));


my $chr_slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome', '6');
my $features = $sfa->fetch_all_by_Slice($chr_slice);
is(@$features, 77, "Found features on chromosome 6");
#print_features($features);


my $ctg_slice = $dba->get_SliceAdaptor->fetch_by_region('contig','AL449423.14');
$features = $sfa->fetch_all_by_Slice($ctg_slice);
debug('-- contig AL449423.14 simple features ---');
is(@$features, 3, "Found 3 features");
print_features($features);

#retrieve a feature via dbID

debug('---- fetch_by_dbID (default coords) ----');
my $feat = $sfa->fetch_by_dbID(49);
is($feat->dbID, 49, "Correct db id");
is($feat->slice->seq_region_name(), '6', "Correct region name");
is($feat->start, 112575158, "Correct feature start");
is($feat->end, 112575740, "Correct feature end");
is($feat->strand, -1, 'Correct feature strand');

print_features([$feat]);

debug('---- transform to clone ----');
$feat = $feat->transform('contig');
is($feat->dbID, 49);
is($feat->slice->seq_region_name(), 'AL590106.7');
is($feat->start, 6339);
is($feat->end, 6921);
is($feat->strand, -1);
print_features([$feat]);


# List_dbidx
my $ids = $sfa->list_dbIDs();
ok (@{$ids});


sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      my $seqname = $f->slice->seq_region_name();
      debug($seqname . ' ' . $f->start().'-'.$f->end().'('.$f->strand().
            ') ['. $f->dbID.'] '. $f->name.' '.$f->gene_id() );
    } else {
      debug('UNDEF');
    }
  }
}

done_testing();

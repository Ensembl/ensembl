use strict;

use lib 't';
use TestUtils qw(test_getter_setter debug);

BEGIN { $| = 1;  
	use Test;
	plan tests => 10;
}

use MultiTestDB;
use Bio::EnsEMBL::SimpleFeature;


our $verbose = 1;

my $multi = MultiTestDB->new;
 
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
print_features($features);

my $cln_slice = $dba->get_SliceAdaptor->fetch_by_region('clone','AL031658.11');
$features = $sfa->fetch_all_by_Slice($cln_slice);

debug('-- cln AL031658.11 simple features ---');
print_features($features);

my $sprctg_slice = $dba->get_SliceAdaptor->fetch_by_region('supercontig',
                                                         'NT_028392');
$features = $sfa->fetch_all_by_Slice($sprctg_slice);

debug('-- sprctg NT_028392 simple features ---');
print_features($features);

my $ctg_slice = $dba->get_SliceAdaptor->fetch_by_region('contig',
                                                       'AL031658.11.1.162976');

$features = $sfa->fetch_all_by_Slice($ctg_slice);

debug('--- contig AL031658.11.1.162976 simple features ---');
print_features($features);


# List_dbidx
my $ids = $sfa->list_dbIDs();
ok (@{$ids});



sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
    my $analysis = $f->analysis->logic_name();
    debug($f->start().'-'.$f->end().'('.$f->strand().') '.
          $f->display_label.' '.$f->score() . " ($analysis)");
  }
}

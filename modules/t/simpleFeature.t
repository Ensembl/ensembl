use strict;

use lib 't';
use TestUtils qw(test_getter_setter);

BEGIN { $| = 1;  
	use Test;
	plan tests => 13;
}

use MultiTestDB;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::RawContig;
use Bio::Seq;


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
# 8 attach a contig
#
# create a dummy seq and contig
#
my $seq  = Bio::Seq->new(-seq => 'ATGCAGCTAGCATCGATGACATCG',
                         -id => 'dummy_contig',
                         -accession => 'dummy_contig');
  
my $contig = Bio::EnsEMBL::RawContig->new();
 
my $name =  'dummy_contig';
$contig->id($name);
$contig->embl_offset(0);
$contig->seq($seq);

# now attach the contig

$sf->contig($contig);
ok($sf);


#
# 9 check adaptor attaching
#
$sf->adaptor($sfa);
ok($sf->adaptor->isa('Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor'));

# list_dbIDs
my $ids = $sfa->list_dbIDs();
ok (@{$ids});


#
# 10 test store method
#
$multi->hide('core', 'simple_feature');


my $analysis_adaptor = $dba->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name('cpg');

my $contig_adaptor = $dba->get_RawContigAdaptor();
$contig = $contig_adaptor->fetch_by_name('AL031658.11.1.162976');

my $simple_feature = Bio::EnsEMBL::SimpleFeature->new
	(-start => 10,
         -end => 20,
         -strand => 1,
         -score  => 20,
         -analysis => $analysis);

$simple_feature->contig($contig);
$simple_feature->display_label('test');


$sfa->store($simple_feature);

my $sth = $dba->prepare('SELECT simple_feature_id from simple_feature');
$sth->execute();
ok($sth->rows() == 1);
my ($id) = $sth->fetchrow_array();

ok($id == $simple_feature->dbID());
ok($simple_feature->adaptor == $sfa);

$multi->restore();


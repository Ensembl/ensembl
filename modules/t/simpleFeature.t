use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 5;
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
$sf = new Bio::EnsEMBL::SimpleFeature;
$sf->start(10);
$sf->end  (14);
$sf->strand(1);
$sf->score(42);
ok($sf);


#
# 2 attach a contig
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
# 3 display_label
#

$sf->display_label('dummy-label');
ok($sf->display_label eq 'dummy-label');


#
# 4 dbID
#
$sf->dbID('42');
ok($sf->dbID == 42);


#
# 5 adaptor
#
$sf->adaptor($sfa);
ok($sf->adaptor->isa('Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor'));
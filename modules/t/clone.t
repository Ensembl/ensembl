use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 8;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use EnsTestDB;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::PerlDB::Clone;
use Bio::EnsEMBL::PerlDB::Contig;



$loaded = 1;

ok(1);

# Database will be dropped when this
# object goes out of scope
my $ens_test = EnsTestDB->new;

$ens_test->do_sql_file("t/minidatabase.dump");

ok($ens_test);



my $db = $ens_test->get_DBSQL_Obj;

$clone_ad = $db->get_CloneAdaptor();


my $clone = Bio::EnsEMBL::PerlDB::Clone->new();

$clone->id('dummy_clone');
$clone->embl_id('dummy_clone');
$clone->embl_version(1);
$clone->version(1);
$clone->htg_phase(3);


ok($clone);

my $seq  = Bio::Seq->new(-seq => 'ATGCAGCTAGCATCGATGACATCG',
		         -id => 'dummy_contig',
	              	 -accession => 'dummy_contig');


ok($seq);



my $contig = Bio::EnsEMBL::PerlDB::Contig->new();


my $name =  'dummy_contig';
$contig->id($name);
$contig->seq($seq);
$contig->order(1);
$contig->embl_offset(0);
$contig->international_name($name);

$clone->add_Contig($contig);

print STDERR "CONTIG ".$contig->primary_seq()."\n";
ok($contig);

print STDERR "contig internaltional name = ".$contig->international_name."\n";

$clone_ad->store($clone);

ok(1);


my $clone_out = $clone_ad->fetch_by_name('dummy_clone');

ok($clone_out);

ok($clone_out->embl_version == 1);

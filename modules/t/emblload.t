
# $Id$
# testing of translations of exons that lie across contig boundaries.
# based on staticgoldenpath.t and staticgoldenpath.dump

## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..4\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::EMBLLOAD::Obj;

use Bio::SeqIO;

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.

$" = ", ";                          # for easier list-printing
    
my $ens_test = EnsTestDB->new();
my $db = $ens_test->get_DBSQL_Obj;



$file = "t/roa1.dat";

$seqio = Bio::SeqIO->new( '-format' => 'EMBL',-file => $file);

while( my $seq = $seqio->next_seq ) {
    
    $obj = Bio::EnsEMBL::EMBLLOAD::Obj->new(-seq => $seq);
    
    ($clone) = $obj->get_Clone();

    $db->write_Clone($clone);

    @genes = $clone->get_all_Genes();

    foreach $gene ( @genes ) {
	$db->write_Gene($gene);
    }

}

print "ok 2\n";

$clone = $db->get_Clone('HSHNRNPA');

print "ok 3\n";

$gene_obj =$db->gene_Obj();

$gene = $gene_obj->get('HSHNRNPA.gene.1');

print "ok 4\n";

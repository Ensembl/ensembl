
# $Id$
# testing of translations of exons that lie across contig boundaries.
# based on staticgoldenpath.t and staticgoldenpath.dump

## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..5\n";
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

my $gadp = $db->get_GeneAdaptor();
my $id;
my %chash;

while( my $seq = $seqio->next_seq ) {
    
    $obj = Bio::EnsEMBL::EMBLLOAD::Obj->new(-seq => $seq);
    
    ($clone) = $obj->get_Clone();

    $db->write_Clone($clone);



    @genes = $clone->get_all_Genes();

    foreach $gene ( @genes ) {
	# make this a subroutine
	foreach $exon ( $gene->get_all_Exons ) {
	     if( !exists $chash{$exon->seqname} ) {
	           $chash{$exon->seqname} = $db->get_Contig($exon->seqname);
             }
             $exon->contig_id($chash{$exon->seqname}->internal_id);
        }
	$id = $gadp->store($gene);
    }

}

print "ok 2\n";

$clone = $db->get_Clone('HSHNRNPA');

print "ok 3\n";


$gene = $gadp->fetch_by_dbID($id);


print "ok 4\n";

if( scalar($gene->get_all_Exons) == 9 ) {
   print "ok 5\n";
} else {
   print "not ok 5\n";
}
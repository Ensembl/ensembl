

use Test;
BEGIN { plan tests => 5 }



use lib './t';
use EnsTestDB;
use Bio::EnsEMBL::DBSQL::GenomicAlignAdaptor;
use Bio::AlignIO;

ok 1;

    
my $ens_test = EnsTestDB->new();

$ens_test->do_sql_file("t/genomicalign.dump");

my $db = $ens_test->get_DBSQL_Obj;

ok $ens_test;


$gadp = Bio::EnsEMBL::DBSQL::GenomicAlignAdaptor->new($db);

ok $gadp;

my $align = $gadp->fetch_GenomicAlign_by_dbID(1);

ok $align;

$alignblockset = $align->get_AlignBlockSet(1);

($a,$alignblock) = $alignblockset->get_AlignBlocks;

ok ($alignblock->start == 25 && $alignblock->end == 30 && $alignblock->align_start == 15 
	&& $alignblock->align_end == 20);

$alignout = Bio::AlignIO->new( -format => 'fasta',-fh => \*STDERR );

$alignout->write_aln($align);



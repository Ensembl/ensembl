#Contact: Emmanuel Mongin (mongin@ebi.ac.uk)

use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::SeqIO;



my $host      = 'ecs2d';
my $dbuser    = 'ecs2dadmin';
my $dbname    = 'caenorhabditis_elegans_core_10_93';
my $dbpass    = 'TyhRv';
my $port;

my $file = "/acari/work4/mongin/mapping/elegans/Primary/wormpep.table";
my %map;

print STDERR "Connecting to $host, $dbname\n";


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-port'   => $port
					    );

open (IN,"$file") || die "can't open file\n";

while (<IN>) {
    chomp;

    my @array = split;

    $map{$array[0]} = $array[1];

}

my $adaptor = $db->get_DBEntryAdaptor();

my $query = "select t.translation_id, ts.stable_id, gs.stable_id from transcript t, gene_stable_id gs, transcript_stable_id ts where t.gene_id = gs.gene_id and t.transcript_id = ts.transcript_id";
my $sth = $db->prepare($query);
$sth->execute();
while (my @res = $sth->fetchrow) {
    my $transl_dbid = $res[0];
    my $transc_stable_id = $res[1];
    my $gene_stable_id = $res[2];

    my $db1 = "wormbase_gene";
    my $db2 = "wormbase_transcript";
    my $db3 = "wormpep_id";
       
    #print  "$gene_stable_id\t$transc_stable_id\n";

    my $dbentry = Bio::EnsEMBL::DBEntry->new
	( -adaptor => $adaptor,
	  -primary_id => $gene_stable_id,
	  -display_id => $gene_stable_id,
	  -version => 1,
	  -release => 1,
	  -dbname => $db1);
    $dbentry->status("XREF");
	
#	print STDERR "$db1\t$db2\n";
    
#    $adaptor->store($dbentry,$transl_dbid,"Translation");
    
    my $transdbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $transc_stable_id,
	      -display_id => $transc_stable_id,
	      -version => 1,
	      -release => 1,
	      -dbname => $db2);
    $transdbentry->status("XREF");
    
#    $adaptor->store($transdbentry,$transl_dbid,"Translation");

    my $ce = $map{$transc_stable_id};
    
    if ($ce) {
	print STDERR "$db3\t$ce\n";
	my $ceentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $ce,
	      -display_id => $ce,
	      -version => 1,
	      -release => 1,
	      -dbname => $db3);
	$ceentry->status("XREF");
	$adaptor->store($ceentry,$transl_dbid,"Translation");
    }

}

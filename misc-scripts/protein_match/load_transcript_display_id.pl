#Contact: Emmanuel Mongin (mongin@ebi.ac.uk)

use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::SeqIO;
#use MultiTestDB;


BEGIN {
    my $script_dir = $0;
    $script_dir =~ s/(\S+\/)\S+/$1/;
    unshift (@INC, $script_dir);
    require "mapping_conf.pl";
}

my %conf =  %::mapping_conf; # configuration options


# global vars


my $dbname     = $conf{'db'};
my $host       = $conf{'host'};
my $user       = $conf{'dbuser'};
my $pass       = $conf{'password'};
my $port       = $conf{'port'};
my $organism   = $conf{'organism'};
my %priority;

$priority{'BRIGGSAE_HYBRID'} = 1000;
$priority{'HUGO'} = 1000;
$priority{'MarkerSymbol'} = 1000;
$priority{'wormbase_transcript'} = 1000;
$priority{'flybase_symbol'} = 1000;
$priority{'Anopheles_symbol'} = 1000;
$priority{'SWISSPROT'} = 900;
$priority{'RefSeq'} = 800;
$priority{'SPTREMBL'} = 700;
$priority{'LocusLink'} = 100;


if (!defined $organism) {    die "\nSome basic options have not been set up, have a look at mapping_conf\nCurrent set up (required options):\norganism: $organism\n\n";
}

print STDERR "Connecting to the database...$dbname $host\n";
print STDERR "dealing with organism ".$organism."\n";


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => $user,
        -dbname => $dbname,
        -host   => $host,
	-pass   => $pass,
	-port   => $port,
        -driver => 'mysql',
	);


my $transadaptor = $db->get_TranscriptAdaptor();
my $geneadaptor  = $db->get_GeneAdaptor();
my $xrefadaptor  = $db->get_DBEntryAdaptor();

my $query = "select transcript_id from transcript";
my $sth = $db->prepare($query);
$sth->execute();

print STDERR "Getting transcript display xref_id\n";
while(my $id = $sth->fetchrow) {
    my $trans = $transadaptor->fetch_by_dbID($id);
    my $xrefs = $trans->get_all_DBLinks;
    my $display;
    my $current = 0;
    my $cdb;
    foreach my $xref(@$xrefs) {
	
	if ($priority{$xref->database} > $current) {
	    $display = $xref;
	    $cdb = $xref->database;
	    $current = $priority{$xref->database};
	}
    }

    $trans->display_xref($display);
    $transadaptor->update($trans);
}

print STDERR "Done\n";
print STDERR "Getting gene display_xref_id\n";


if ($organism ne "elegans") {
    my $query1 = "select gene_id from gene";
    my $sth1 = $db->prepare($query1);
    $sth1->execute();
    while(my $gene_id = $sth1->fetchrow) {
	my $gene = $geneadaptor->fetch_by_dbID($gene_id);
	my $transcripts = $gene->get_all_Transcripts();
	my $display;
	my $current;
	TRANS:foreach my $trans(@$transcripts) {
	    my $xref = $transadaptor->get_display_xref($trans);
	    if(!$xref){
	      next TRANS;
	    }
	    if ($priority{$xref->database} > $current) {
	      $display = $xref;
	      $current = $priority{$xref->database};
	    }
	}
      
	$gene->display_xref($display);
	$geneadaptor->update($gene);
       
    }
}
#Not sure id it is really needed if wormbase_gene is put in the priority list... Laura?
elsif ($organism eq "elegans") {
    my $query1 = "select g.gene_id, x.xref_id from gene_stable_id g, xref x, external_db e where g.stable_id = x.display_label and x.external_db_id = e.external_db_id and e.db_name = 'wormbase_gene'";
    my $sth1 = $db->prepare($query1);
    $sth1->execute();
    
    while(my ($gene_id,$xref)  = $sth1->fetchrow) {
      my $xref_obj = $db->get_DBEntryAdaptor->fetch_by_dbID($xref);
      my $gene = $geneadaptor->fetch_by_dbID($gene_id);
	$gene->display_xref($xref_obj);
	$geneadaptor->update($gene);
     }
    my $query1 = "select g.gene_id, x.xref_id from gene_stable_id g, xref x, external_db e where g.stable_id = x.display_label and x.external_db_id = e.external_db_id and e.db_name = 'wormbase_pseudogene'";
    my $sth1 = $db->prepare($query1);
    $sth1->execute();
    
    while(my ($gene_id,$xref)  = $sth1->fetchrow) {
      my $xref_obj = $db->get_DBEntryAdaptor->fetch_by_dbID($xref);
      my $gene = $geneadaptor->fetch_by_dbID($gene_id);
	$gene->display_xref($xref_obj);
	$geneadaptor->update($gene);
     }

     my $query1 = "select t.transcript_id, x.xref_id from transcript_stable_id t, xref x, external_db e where t.stable_id = x.display_label and x.external_db_id = e.external_db_id and e.db_name = 'wormbase_pseudogene'";
    my $sth1 = $db->prepare($query1);
    $sth1->execute();
    
    while(my ($gene_id,$xref)  = $sth1->fetchrow) {
      my $xref_obj = $db->get_DBEntryAdaptor->fetch_by_dbID($xref);
      my $trans = $transadaptor->fetch_by_dbID($gene_id);
	$trans->display_xref($xref_obj);
	$transadaptor->update($trans);
     }
    
}

print STDERR "Done\n";




#Contact: Emmanuel Mongin (mongin@ebi.ac.uk)
#This script loads the display_xref ids in the tables gene and transcript. The priority to which xref to load is given in the priority hash. Please edit it if needed. These script can be used for any organism
#The database connection values are taken from mapping_conf.pl



use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::SeqIO;

BEGIN {
    my $script_dir = $0;
    $script_dir =~ s/(\S+\/)\S+/$1/;
    unshift (@INC, $script_dir);
    require "mapping_conf.pl";
}

my %conf =  %::mapping_conf; # configuration options


# Database connection values
my $dbname     = $conf{'db'};
my $host       = $conf{'host'};
my $user       = $conf{'dbuser'};
my $pass       = $conf{'password'};
my $organism   = $conf{'organism'};
my %priority;

#Priority set up
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
$priority{'Anopheles_paper'} = 50;
$priority{'Celera_Gene'} = 50;

if (!defined $organism) {
    die "\nSome basic options have not been set up, have a look at mapping_conf\nCurrent set up (required options):\norganism: $organism\n\n";
}

print STDERR "Connecting to the database...\n";

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => $user,
        -dbname => $dbname,
        -host   => $host,
	-pass   => $pass,			     
        -driver => 'mysql',
	);

my $transadaptor = $db->get_TranscriptAdaptor();
my $geneadaptor  = $db->get_GeneAdaptor();
my $xrefadaptor  = $db->get_DBEntryAdaptor();

#First get through all of the transcripts and set up the display xref following the priority list
my $query = "select transcript_id from transcript";
my $sth = $db->prepare($query);
$sth->execute();

while(my $id = $sth->fetchrow) {
    my $trans = $transadaptor->fetch_by_dbID($id);
    my $xrefs = $trans->get_all_DBLinks;
    my $display;
    my $current = 0;
    foreach my $xref(@$xrefs) {
	if ($priority{$xref->database} > $current) {
	    $display = $xref->dbID;
	}
    }
    $trans->display_xref($display);
    $transadaptor->update($trans);

#Then get through all of the genes and choose from the transcripts which one to use. There is currently am exeption for elegans as the display xref should be the gene and transcripts stable ids.
if ($organism ne "elegans") {

    my $query1 = "select gene_id from gene";
    my $sth1 = $db->prepare($query1);
    $sth1->execute();
    
    while(my $gene_id = $sth1->fetchrow) {
	my $gene = $geneadaptor->fetch_by_dbID($gene_id);
	my $transcripts = $gene->get_all_Transcripts();
	my $display;
	my $current;
	foreach my $trans(@$transcripts) {
	    my $id = $trans->dbID();
	    my $xrefid = $transadaptor->get_display_xref_id($id);
	    
	    eval {
		my $xref = $xrefadaptor->fetch_by_dbID($xrefid);
		if ($priority{$xref->database} > $current) {
		    $display = $xref->dbID;
		}
	    };
	}
	$gene->display_xref($display);
	$geneadaptor->update($gene);
    }
}
elsif ($organism eq "elegans") {
    my $query1 = "select g.gene_id, x.xref_id from gene_stable_id g, xref x, external_db e where g.stable_id = x.display_label and x.external_db_id = e.external_db_id and e.db_name = 'wormbase_gene'";
    my $sth1 = $db->prepare($query1);
    $sth1->execute();
    
    while(my ($gene_id,$xref)  = $sth1->fetchrow) {
	my $gene = $geneadaptor->fetch_by_dbID($gene_id);
	$gene->display_xref($xref);
	$geneadaptor->update($gene);
     }
    
}






#Contact: Emmanuel Mongin (mongin@ebi.ac.uk)
#This script is supposed to retrieve basic stats about the mapping. More stuff need to be added but at least that's a beginning...

use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

BEGIN {
    my $script_dir = $0;
    $script_dir =~ s/(\S+\/)\S+/$1/;
    unshift (@INC, $script_dir);
    require "mapping_conf.pl";
}

my %conf =  %::mapping_conf; # configuration options


# global vars

#Global variables from the mapping configuration file
my $dbname     = $conf{'db'};
my $host       = $conf{'host'};
my $user       = $conf{'dbuser'};
my $pass       = $conf{'password'};
my $organism   = $conf{'organism'};
my $file       = $conf{'statistic_file'};
my $t_idt      = $conf{'min_ensembl_idt'};
my $q_idt      = $conf{'min_known_idt'};
#print STDERR "have query identity ".$q_idt."\n";
my $unique = undef;
my $refid  = undef;

#Hashes definition
my %dbid2name;

#open (STDOUT,'>$file') or die "couldn't open ".$file;
#print STDERR "have opened ".$file."\n";
#Get the DB object
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => $user,
        -dbname => $dbname,
        -host   => $host,
	-pass   => $pass,			     
        -driver => 'mysql',
	);


#First get the DBid for each dbname
my $query = "select external_db_id, db_name from external_db";

my $sth = $db->prepare($query);

$sth->execute;

while (my $hash = $sth->fetchrow_hashref()) {
    my $dbid = $hash->{'external_db_id'};
    my $dbname = $hash->{'db_name'};
    
    $dbid2name{$dbid} = $dbname;
    if ($dbname eq "RefSeq") {
	$refid = $dbid;
    }
}

my @ids = keys %dbid2name;

foreach my $i(@ids) {
    $query = "select count(distinct(dbprimary_acc)) from xref where external_db_id = $i";
    $sth = $db->prepare($query);
    $sth->execute;
    $unique = $sth->fetchrow;

    print STDOUT "$unique unique $dbid2name{$i} entries matched\n" if($unique);
}



$query = "select count(distinct(ensembl_id)) from object_xref";

$sth = $db->prepare($query);
$sth->execute;
$unique = $sth->fetchrow;

print STDOUT "$unique Ensembl id matched\n";

print STDOUT "REFID: $refid\n";

while ($q_idt < 100) {
    
    $query = "select count(distinct(o.object_xref_id)) 
          from object_xref as o,
          xref as x,
          identity_xref as i
          where  o.xref_id = x.xref_id
          and external_db_id = $refid
          and o.object_xref_id = i.object_xref_id
          and i.query_identity >= $q_idt";
    #print $query."\n";
    $sth = $db->prepare($query);
    $sth->execute;
    my $row = $sth->fetchrow();
    print STDOUT "$q_idt\t$row\n";
    $q_idt = $q_idt + 5;
    if ($q_idt > 100) {
	$q_idt = 100;
    }
}
$sth->finish;
#close(OUT);

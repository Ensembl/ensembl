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
my $t_idt      = $conf{'target_idt'};
my $q_idt      = $conf{'query_idt'};

my $unique = undef;
my $refid  = undef;

#Hashes definition
my %dbid2name;

open (OUT,'>$file');

#Get the DB object
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => $user,
        -dbname => $dbname,
        -host   => $host,
	-pass   => $pass,			     
        -driver => 'mysql',
	);


#First get the DBid for each dbname
my $query = "select externalDBId, db_name from externalDB";

my $sth = $db->prepare($query);

$sth->execute;

while (my $hash = $sth->fetchrow_hashref()) {
    my $dbid = $hash->{'externalDBId'};
    my $dbname = $hash->{'db_name'};
    
    $dbid2name{$dbid} = $dbname;
    if ($dbname eq "RefSeq") {
	$refid = $dbid;
    }
}

my @ids = keys %dbid2name;

foreach my $i(@ids) {
    $query = "select count(distinct(dbprimary_id)) from Xref where externalDBId = $i";
    $sth = $db->prepare($query);
    $sth->execute;
    $unique = $sth->fetchrow;

    print OUT "$unique unique $dbid2name{$i} entries matched\n";
}



$query = "select count(distinct(ensembl_id)) from objectXref";

$sth = $db->prepare($query);
$sth->execute;
$unique = $sth->fetchrow;

print OUT "$unique Ensembl id matched\n";

print OUT "REFID: $refid\n";

while ($q_idt < 100) {
    
    $query = "select count(distinct(o.objectxrefId)) 
          from objectXref as o,
          Xref as x,
          identityXref as i
          where  o.xrefId = x.xrefId
          and externalDBId = $refid
          and o.objectxrefId = i.objectxrefId
          and i.query_identity >= $q_idt";

    $sth = $db->prepare($query);
    $sth->execute;
    my $row = $sth->fetchrow();
    print OUT "$q_idt\t$row\n";
    $q_idt = $q_idt + 5;
    if ($q_idt > 100) {
	$q_idt = 100;
    }
}

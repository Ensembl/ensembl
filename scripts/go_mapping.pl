use strict;
use Bio::EnsEMBL::DBLoader;
use GD;

$| = 1;

my $iprgo = shift(@ARGV);
open (IPRGO,"<$iprgo");
while (<IPRGO>) {
    my ($ipr,$go) = split (/\t/);
    print STDERR "INTERPRO $ipr maps to GO $go\n";
}
exit;
#DB parameters
my $dbtype = 'rdb';
my $host   = 'ecs1c';
my $port   = '';
my $dbname = 'ensembl080';
my $dbuser = 'ensadmin';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";

print STDERR "Using $locator for db\n";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);
my $sp = $db->get_StaticGoldenPathAdaptor;

my $sth=$db->prepare("select distinct interpro_ac from interpro");
$sth->execute();
print STDERR "Getting all interpro domains...\n";
my @domains;
while (my $row= $sth->fetchrow_array) {
#    print STDERR "Getting $row\n";
    push (@domains,$row);
}

foreach my $domain (@domains) {
    print STDERR "Processing domain $domain\n";
    print STDERR "Finding genes...\n";
    my @genes = find_Genes($db,$domain);
    
    print STDERR "Processing genes...\n";
    my %bychr=process_Genes($domain,@genes);

    print STDERR "Dumping genes...\n";
    dump_genes($domain,@genes);
}

sub find_Genes {
    my ($db,$domain) = @_;

    my @genes;
    
    #Find all genes that have an interpro domain, 
    #and count how many domains it has in total
    my $sth = $db->prepare("select tr.gene from interpro i,protein_feature pf,transcript tr where i.interpro_ac = '$domain' and i.id = pf.hid and pf.translation = tr.translation group by tr.gene");

    $sth->execute();
    my $exnum;
    while (my $rowhash = $sth->fetchrow_hashref) {
	my $geneid = $rowhash->{'gene'};
	
    }
    return @genes;
}

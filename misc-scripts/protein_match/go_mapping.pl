use strict;
use Bio::EnsEMBL::DBLoader;
use GD;

$| = 1;

my %go_map;
my $iprgo = shift(@ARGV);
open (IPRGO,"<$iprgo");
while (<IPRGO>) {
    chomp();
    my ($ipr,$go) = $_ =~ /(IPR\S+)\|GO\:(\d+)/;
    $go_map{$ipr}=$go;
}

#DB parameters
my $dbtype = 'rdb';
my $host   = 'ecs1b';
my $port   = '';
my $dbname = 'ensembl100';
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
print STDERR "Getting genes for each domain\n";
foreach my $domain (@domains) {
    my @genes;
    
    #Find all genes that have an interpro domain, 
    #and count how many domains it has in total
    my $sth = $db->prepare("select tr.gene from interpro i,protein_feature pf,transcript tr where i.interpro_ac = '$domain' and i.id = pf.hid and pf.translation = tr.translation group by tr.gene");

    $sth->execute();
    my $exnum;
    while (my $rowhash = $sth->fetchrow_hashref) {
	my $geneid = $rowhash->{'gene'};
	if ($go_map{$domain}) {
	    print "$geneid\t".$go_map{$domain}."\n";
	}
    }
}







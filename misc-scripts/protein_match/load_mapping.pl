use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;

my %hugosyn;
my %hugosymbol;
my %scopsyn;
my %gene_map;
my %transcript_map;
my ($mapping, $hugosyn, $scopsyn, $out);



&GetOptions(
	    'mapping:s'=>\$mapping,
	    'hugosyn:s'=>\$hugosyn,
	    'scopsyn:s'=>\$scopsyn
            );

my $dsn = "DBI:mysql:database=ensembl090_tmp;host=ecs1c";    
my $db = DBI->connect("$dsn",'ensadmin') || die ("Could not connect to db!");   

my $adaptor =  Bio::EnsEMBL::DBSQL::DBEntryAdaptor->new($db);


#open (MAPS, "$map");
#while (<MAP>) {
#    chomp;
#    my ($transcript,$gen
#Read Hugo file to get out synonyms

open (HUGO, "$hugosyn") || die "Can't open file $mapping\n";
while (<HUGO>) {
    chomp;
    my ($hgnc, $symbol, $alias, $withdrawn) = split (/\t/,$_);
    
    my @aliases = split (/, /,$alias);
    my @withdrawns = split (/, /,$withdrawn);
    
    $hugosymbol{$symbol}=$hgnc;

    foreach my $al(@aliases) {
	push(@{$hugosyn{$symbol}},$al);
    }

    foreach my $wi(@withdrawns) {
	push(@{$hugosyn{$symbol}},$wi);
    }
}
close (HUGO);


#Read SCOP file to get out synonyms
open (SCOP, "$scopsyn") || die "Can't open file $scopsyn\n";
while (<SCOP>) {
    chomp;
    my ($scopac, $pdb, $chain, $scopnb) = split(/\t/,$_);
    
    #my $uni = "$pdb||$chain";

    push(@{$scopsyn{$scopac}},$pdb);
    push(@{$scopsyn{$scopac}},$chain);
    push(@{$scopsyn{$scopac}},$scopnb);
}
close (SCOP);

#Read final mapping
open (MAPPING, "$mapping") || die "Can't open file $mapping\n";
while (<MAPPING>) {
    chomp;
    my ($ens, $db, $primary_ac) = split(/\t/,$_);

#Get SP mapping
    if (($db ne "HUGOSYMBOL") && ($db ne "SCOP") && ($db ne "SCOP1") && ($db ne "HUGOID") && ($db ne "HUGOALIAS") && ($db ne "HUGOWITHDRAWN")) {
	my ($ac1) = $ens =~ /COBP(\d+)/;
	$ens = "COBT"."$ac1";
	my $dbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $primary_ac,
	      -display_id => $primary_ac,
	      -version => 1,
	      -release => 1,
	      -dbname => $db );
	$adaptor->store($dbentry,$ens,"Gene");
    }

    if ($db eq "HUGOSYMBOL") {
	#print STDERR "HERE\n";
	my ($ac1) = $ens =~ /COBP(\d+)/;
	$ens = "COBT"."$ac1";
	
	my $dbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $primary_ac,
	      -display_id => $primary_ac,
	      -version => 1,
	      -release => 1,
	      -dbname => $db );
	if ($hugosyn{$primary_ac}) {
	    my @synonyms = @{$hugosyn{$primary_ac}};
	    #print STDERR "SYN: @synonyms\n";
	    foreach my $syn (@synonyms) {
		
		if ($syn =~ /\S+/) {
		    #print STDERR "$syn\n";
		    $dbentry->add_synonym($syn);
		   
		}
	    }
	}
	
	
	$adaptor->store($dbentry,$ens,"Gene");
    }

    if ($db eq "SCOP") {
     	my ($ac1) = $ens =~ /COBP(\d+)/;
	$ens = "COBT"."$ac1";
	my $dbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $primary_ac,
	      -display_id => $primary_ac,
	      -version => 1,
	      -release => 1,
	      -dbname => $db );
	if ($scopsyn{$primary_ac}) {
	my @synonyms = @{$scopsyn{$primary_ac}};
	foreach my $syn (@synonyms) {
	    if ($syn =~ /\S+/) {
			$dbentry->add_synonym($syn);
	    }
	}
    }
	$adaptor->store($dbentry,$ens,"Gene");
	
    }

}

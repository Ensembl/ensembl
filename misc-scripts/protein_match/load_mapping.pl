use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::SeqIO;

my %hugosyn;
my %hugoid;
my %scopsyn;
my %scopid;
my %gene_map;
my %transcript_map;
my %spid;
my ($mapping, $hugosyn, $scopsyn, $out, $spsyn);

#perl ../../src/ensembl-live/misc-scripts/protein_match/load_mapping.pl -mapping outputs/final_sorted.map -hugosyn secondary/ens4.txt -scopsyn secondary/dir.dom.scop.txt_1.53 -spid primary/hum_sp_sptrembl.pep

&GetOptions(
	    'mapping:s'=>\$mapping,
	    'hugosyn:s'=>\$hugosyn,
	    'scopsyn:s'=>\$scopsyn,
	    'spid:s'=>\$spsyn
            );

my $dsn = "DBI:mysql:database=ensembl090_tmp;host=ecs1c";    
my $db = DBI->connect("$dsn",'ensadmin') || die ("Could not connect to db!");   

my $adaptor = Bio::EnsEMBL::DBSQL::DBEntryAdaptor->new($db);



print STDERR "Getting SP mapping\n";

#my $in  = Bio::SeqIO->new(-file => $spsyn, '-format' =>'swiss');

#while ( my $seq = $in->next_seq() ) {
    my $ac;# = $seq->accession;
    my $id;# = $seq->id;
    $spid{$ac} = $id;
#}


open (HUGO, "$hugosyn") || die "Can't open file $mapping\n";
while (<HUGO>) {
    chomp;

#get red of the cariage return present in Hugos
    $_ =~ s/\r//g;
    my ($hgnc, $symbol, $alias, $withdrawn) = split (/\t/,$_);
    
    my @aliases = split (/, /,$alias);
    my @withdrawns = split (/, /,$withdrawn);
    
    $hugoid{$hgnc}=$symbol;

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

#Set up the display id
    my $display = $pdb." ".$chain;

    push (@{$scopid{$scopac}},$display);

    #push(@{$scopsyn{$scopac}},$pdb);
    #push(@{$scopsyn{$scopac}},$chain);

#Scop number becomes a synonym (not stable)
    push(@{$scopsyn{$scopac}},$scopnb);
}
close (SCOP);

#Read final mapping
open (MAPPING, "$mapping") || die "Can't open file $mapping\n";
while (<MAPPING>) {
    chomp;
    $_ =~ s/\r//g;
    my ($ens, $db, $primary_ac) = split(/\t/,$_);
    
#Get SP mapping
    #if (($db ne "HUGOSYMBOL") && ($db ne "SCOP") && ($db ne "SCOP1") && ($db ne "SCOP2") && ($db ne "HUGOID") && ($db ne "HUGOALIAS") && ($db ne "HUGOWITHDRAWN")) {
    if (($db eq "EMBL") || ($db eq "EC") || ($db eq "OMIM") || ($db eq "REFSEQ") || ($db eq "LOCUS")) {
	

##############Temporary changes###########################
	my ($ac1) = $ens =~ /COBP(\d+)/;
	$ens = "COBT"."$ac1";
##########################################################	
	
	
	my $dbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $primary_ac,
	      -display_id => $primary_ac,
	      -version => 1,
	      -release => 1,
	      -dbname => $db );
	$adaptor->store($dbentry,$ens,"Gene");
    }
    
    

    if (($db eq "SP") || ($db eq "SPTREMBL")) {

	if (!defined $spid{$primary_ac}) {
	    #print "SP primary Ac ($primary_ac) does not have an id\n";
	} 
	
	my $dbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $primary_ac,
	      -display_id => $spid{$primary_ac},
	      -version => 1,
	      -release => 1,
	      -dbname => $db );
    }
    
    if ($db eq "HUGOID") {

##################Temporary changes#######################
	my ($ac1) = $ens =~ /COBP(\d+)/;
	$ens = "COBT"."$ac1";
##########################################################
	
	if (!defined $hugoid{$primary_ac}) {
	    print "Hugo primary Ac ($primary_ac) does not have an id\n";
	} 


	my $dbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $primary_ac,
	      -display_id => $hugoid{$primary_ac},
	      -version => 1,
	      -release => 1,
	      -dbname => $db );
	if ($hugosyn{$primary_ac}) {
	    my @synonyms = @{$hugosyn{$primary_ac}};
	    foreach my $syn (@synonyms) {
		
		if ($syn =~ /\S+/) {
		    $dbentry->add_synonym($syn);
		}
	    }
	}
	
	
	$adaptor->store($dbentry,$ens,"Gene");
    }

if ($db eq "SCOP") {
     	
#############tmp########################
    my ($ac1) = $ens =~ /COBP(\d+)/;
    $ens = "COBT"."$ac1";
########################################

    if (!defined $scopid{$primary_ac}) {
	    print "SCOP primary Ac ($primary_ac) does not have an id\n";
	} 

	my $dbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $primary_ac,
	      -display_id => $scopid{$primary_ac},
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




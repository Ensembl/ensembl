use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $host      = 'ecs1b';
my $dbuser    = 'ensro';
my $dbname    = 'anopheles_arne_core_9_2';
my $dbpass    = undef;
my $path      = 'NOTRE_DAME';
my $port;
my $mapping;
my $matepairs;
my %map;
my %seen;

#perl format_matepairs.pl -mapping /acari/work4/mongin/anopheles_gambiae/mapacs/celera2uid2embl.txt -matepairs /acari/work4/mongin/anopheles_gambiae/matepairs/mosquito_scaffold_matepair.txt


&GetOptions(
	    'db:s' => \$dbname,
	    'dbhost:s'=> \$host,
	    'dbuser:s'=> \$dbuser,
	    'mapping:s'=>\$mapping,
	    'matepairs:s'=>\$matepairs
	    );

print STDERR "Connecting to $host, $dbname\n";

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-port'   => $port
					   );



open (MAP,"$mapping") || die "Can't open mapping (celera uid to celera stable ids): $mapping\n";

while (<MAP>) {
    chomp;
    my @array = split;
    
    $map{$array[1]} = $array[2];

}

close (MAP);

open (MATE,"$matepairs") || die "Can't open matepairs file: $matepairs\n";

while (<MATE>) {
    chomp;
    my @array = split;
    my $cra_id;
    
    if ($map{$array[0]}) {
	$cra_id = $map{$array[0]};
    }
    else {
	#print STDERR "Celera internal id $array[0] not defined in mapping\n";
    }

    #if ($array[12] ne "happy") {
	my $tag = $array[12];


	#There is some cases where primer start = primer end, take these cases away
	if ($array[3] == $array[4]) {
	    next;
	}
	
	if ($array[7] == $array[8]) {
	    next;
	}

	my @coord = ($array[3],$array[4],$array[7],$array[8]);


	@coord =  sort {$a <=> $b} @coord;
	
	my $start = $coord[0];
	my $end = $coord[3];

	if (($start == $array[4]) || ($start == $array[8])) {
	    die "oups...check the code\n";
	}
	
	
	my ($chr_start,$chr_end,$ori,$chr_id) = &fpc2chr($cra_id,$start,$end);

#	my $chr_id;

	#if (! defined $seen{$cra_id}) {
	#    ($chr_id) = &contig2chr($cra_id);
	#    $seen{$cra_id} = $chr_id;
	#}
	
	print "$chr_id\t$chr_start\t$chr_end\t1\t$tag\n";
    
    #}
}    


sub fpc2chr {
    my ($sc,$start,$end) = @_;
    
    $start = $start + 1;
    $end = $end;

     my $query = "SELECT  
   if(a.superctg_ori=1,($start-a.superctg_start+a.chr_start),
                    (a.chr_start+a.superctg_end-$end)),
   if(a.superctg_ori=1,($end-a.superctg_start+a.chr_start),
                    (a.chr_start+a.superctg_end-$start)),
     c.name, a.superctg_ori FROM assembly a, chromosome c WHERE a.superctg_name = '$sc' and c.chromosome_id = a.chromosome_id";

    my $sth = $db->prepare($query);
    $sth->execute;
    
    my ($chr_start,$chr_end,$name,$ori) = $sth->fetchrow;

    if ($chr_start == 0) {
	print "ID: $sc\tSTART: $start\tEND: $end\tCHR_START: $chr_start\tCHR_END: $chr_end\n";
	print "$query\n";
    }

    return($chr_start,$chr_end,$ori,$name);

}


sub contig2chr {
    my ($cra_id) = @_;

#    my $command = "select ch.name from contig c, chromosome ch where c.contig_id like '$cra_id%' and c.chromosome_id = ch.chromosome_id";

    my $command = "select c.name from chromosome c, assembly a where a.superctg_name = '$cra_id' and a.chromosome_id = c.chromosome_id";

    print STDERR "$command\n";

    my $sth = $db->prepare($command);
    $sth->execute;

    my ($chr_id) = $sth->fetchrow;

    return($chr_id);
}


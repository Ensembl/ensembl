use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $dbname;
my $dbhost;
my $dbuser;
my $dbpass;
my $mapping;
my $matepairs;

&GetOptions(
	    'db:s' => \$dbname,
	    'dbhost:s'=> \$dbhost,
	    'dbuser:s'=> \$dbuser,
	    'mapping:s'=>\$mapping,
	    'matepairs:s'=>\$matepairs
	    );


#perl format_matepairs.pl -db anopheles_gambiae_core_6_1 -dbhost ecs1d -dbuser ensro -mapping /acari/work4/mongin/anopheles_gambiae/matepairs/ga_name_scf_uid_lookup.txt -matepairs /acari/work4/mongin/anopheles_gambiae/matepairs/mosquito_scaffold_matepair.txt

my %map;
my %seen;

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => $dbuser,
        -dbname => $dbname,
        -host   => $dbhost
        );


open (MAP,"$mapping") || die "Can't open mapping (celera uid to celera stable ids): $mapping\n";

while (<MAP>) {
    chomp;
    my @array = split;
    
    $map{$array[1]} = $array[0];

}

close (MAP);

open (MATE,"$matepairs") || die "Can't open matepairs file: $matepairs\n";

while (<MATE>) {
    chomp;
    my @array = split;
    my $cra_id;
    
    my $clone_uid = $array[1];

    if ($map{$array[0]}) {
	$cra_id = $map{$array[0]};
    }
    else {
	#print STDERR "Celera internal id $array[0] not defined in mapping\n";
    }

    if ($array[12] eq "happy") {
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
	
	
	my ($chr_start,$chr_end,$ori) = &fpc2chr($cra_id,$start,$end);

	my $chr_id;

	if (! defined $seen{$cra_id}) {
	    ($chr_id) = &contig2chr($cra_id);
	    $seen{$cra_id} = $chr_id;
	}
	
	my $diff = $chr_end - $chr_start;
	
	if ($diff > 50000) {

	    print "$seen{$cra_id}\t$chr_start\t$chr_end\t1\t$tag\t$clone_uid\n";
	}
    
    }

    

}    


sub fpc2chr {
    my ($cra_id,$start,$end) = @_;
    
    $start = $start + 1;
    $end = $end;

    my $command = "select fpcctg_ori,IF(fpcctg_ori=1,($start+chr_start-fpcctg_start), (chr_start+fpcctg_end-$end)) as chrom_start, IF(fpcctg_ori=1,($end+chr_start-fpcctg_start),  (chr_start+fpcctg_end-$start)) as chrom_end from static_golden_path where fpcctg_name = '$cra_id' and fpcctg_start <= $start and fpcctg_end >= $start";

    my $sth = $db->prepare($command);
    $sth->execute;
    
    my ($ori,$chr_start,$chr_end) = $sth->fetchrow;

    if ($chr_start == 0) {
	print "ID: $cra_id\tSTART: $start\tEND: $end\tCHR_START: $chr_start\tCHR_END: $chr_end\n";
	print "$command\n";
    }
    return($chr_start,$chr_end,$ori);

}


sub contig2chr {
    my ($cra_id) = @_;

    my $command = "select ch.name from contig c, chromosome ch where c.id like '$cra_id%' and c.chromosomeId = ch.chromosome_id";

    my $sth = $db->prepare($command);
    $sth->execute;

    my ($chr_id) = $sth->fetchrow;

    return($chr_id);
}


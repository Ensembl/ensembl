use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBLoader;

my $host = "ecs2d";
my $dbname = "anopheles_gambiae_core_6_1a";

my $enslocator = "Bio::EnsEMBL::DBSQL::DBAdaptor/host=$host;dbname=$dbname;user=ensro;perlonlyfeatures=1";

my $db =  Bio::EnsEMBL::DBLoader->new($enslocator);


my %hash;

open (ID,"/acari/work4/mongin/anopheles_gambiae/mapacs/celera2uid2embl.txt") || die;

while (<ID>) {
    chomp;
    my ($cer,$uid,$embl) = split;

    $hash{$embl} = $cer;

}

open (IN,"/acari/work7a/mongin/assembly/Agambiae_golden_path_17_7_02_XL.txt") || die;
my $current_chr;
my $chr_start;
my $chr_end;
my %seen;
while (<IN>) {
    chomp;
    
    my ($chr,$embl,$start,$end,$ori) = split;
        
    if ($current_chr ne $chr) {
	$chr_start = 0;
	$chr_end = 0;
	$current_chr = $chr;
    }

    foreach my $k (keys %hash) {
	if ($k =~ /$embl$/) {
	    $embl = $hash{$k};
	}
    }
    
    if ($embl) {
	
	my $query = "select MIN(fpcctg_start), MAX(fpcctg_end) from static_golden_path where fpcctg_name = '$embl'";
	my $sth = $db->prepare($query);
	$sth->execute;
	my ($fstart,$fend) = $sth->fetchrow;
	
	if (! defined $seen{$chr}) {
	    $chr_start = $fstart;
	    $chr_end = $fend;
	    $seen{$chr} = 1;
	}
	else {
	    $chr_start = $chr_end + 10000;
	    $chr_end = $chr_start + $fend - 1;
	}

	if ($ori == -1) {
	    $ori = "-";
	}
	else {
	    $ori = "+";
	}

	if ($embl) {
	    print "$embl\t$chr\t$chr_start\t$chr_end\t$ori\n";
	}
    }
}

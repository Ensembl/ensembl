use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::ProteinAdaptor;

use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

my $host      = 'ecs1b';
my $dbuser    = 'ensro';
my $dbname    = 'anopheles_arne_core_9_2';
my $dbpass    = undef;
my $path      = 'NOTRE_DAME';
my $port;#      = 19322;

print STDERR "Connecting to $host, $dbname\n";


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-port'   => $port
					   );

my $bacs = "/acari/work4/mongin/final_build/bacs/BACS.txt";
my $map = "/acari/work4/mongin/final_build/bacs/bacs2bacs_ends.txt";
my %map;

open (BAC,$bacs) || die "Can't open $bacs\n";
open (MAP,$map) || die "Can't open $map\n";

while (<MAP>) {
    chomp;
    my ($ac,$e1,$e2) = split;
    my $end = join(":",$e1,$e2);
    $map{$ac} = $end;
}

while(<BAC>) {
    chomp;
    my (@array) = split;
    my $ac = $array[0];
    my $sc = $array[2];
    my $start = $array[8];
    my $end = $array[9];
    my $status = $array[11];

#    print STDERR "BO:$sc\t$array[4]\t$start\t$end\n";
   
    my $query = "SELECT  
   if(a.superctg_ori=1,($start-a.superctg_start+a.chr_start),
                    (a.chr_start+a.superctg_end-$end)),
   if(a.superctg_ori=1,($end-a.superctg_start+a.chr_start),
                    (a.chr_start+a.superctg_end-$start)),
     c.name FROM assembly a, chromosome c WHERE a.superctg_name = '$sc' and c.chromosome_id = a.chromosome_id";
#     print STDERR "$query\n";
    my $sth = $db->prepare($query);
    $sth->execute();
    my @res = $sth->fetchrow_array();
    my $end = $map{$ac};
    my ($e1,$e2) = split(/:/,$end);
    print "$res[2]\t$res[0]\t$res[1]\t1\t$status\t$ac\t$e1\t$e2\n";
    
}

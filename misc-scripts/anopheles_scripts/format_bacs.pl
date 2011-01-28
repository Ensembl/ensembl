=head1 NAME

  format_bacs.pl

=head1 SYNOPSIS
 
  format_bacs.pl -db anopheles_arne_core_9_2 -dbhost ecs1b -dbuser ensro -bacs /acari/work4/mongin/final_build/bacs/BACS.txt -map /acari/work4/mongin/final_build/bacs/bacs2bacs_ends.txt

=head1 DESCRIPTION

  Takes a BAC mapping in the following format:
  
180L9   219000001128745 AAAB01008847    19611853        X       3715079 1A      5A      45775   151560  IN      unmapped
150B16  219000001128745 AAAB01008847    19611853        X       3715079 1A      5A      47855   170796  IN      unmapped

   and remap the coordinates in chromosome coordinates in a format ready to go to mapfrag

=head1 CONTACTS

  dev@ensembl.org
  mongin@ebi.ac.uk

=cut

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $host      = '';
my $dbuser    = '';
my $dbname    = '';
my $dbpass    = '';
my $path      = '';
my $port;

my $bacs;
my $map;

&GetOptions(
	    'db:s' => \$dbname,
	    'dbhost:s'=> \$dbhost,
	    'pass:s' => \$dbpass,
            'port:s' => \$dbport,
            'dbuser:s'=> \$dbuser,
	    'mapping:s'=>\$map,
	    'bacs:s'=>\$bacs
	    );

print STDERR "Connecting to $host, $dbname\n";


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-port'   => $port
					   );


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
   
    my $query = "SELECT  
   if(a.superctg_ori=1,($start-a.superctg_start+a.chr_start),
                    (a.chr_start+a.superctg_end-$end)),
   if(a.superctg_ori=1,($end-a.superctg_start+a.chr_start),
                    (a.chr_start+a.superctg_end-$start)),
     c.name FROM assembly a, chromosome c WHERE a.superctg_name = '$sc' and c.chromosome_id = a.chromosome_id";

    my $sth = $db->prepare($query);
    $sth->execute();
    my @res = $sth->fetchrow_array();
    my $end = $map{$ac};
    my ($e1,$e2) = split(/:/,$end);
    print "$res[2]\t$res[0]\t$res[1]\t1\t$status\t$ac\t$e1\t$e2\n";
    
}

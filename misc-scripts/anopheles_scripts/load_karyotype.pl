#!/usr/local/ensembl/bin/perl
=head1 NAME


=head1 SYNOPSIS
 
  load_karyotype.pl -db anopheles_arne_core_9_2 -dbhost ecs1b -dbuser ensro -file /acari/work4/mongin/final_build/assembly/Agambiae_golden_path_19_7_02_XL.txt

=head1 DESCRIPTION
 
  This script takes a flat file containing the assembly of the scaffolds as well as the corresponding bands for each scaffolds. The format should be the following:
  chromosome\tscaffold\tstart_band\tend_band\tori
  This loads the karyotype table

=head1 CONTACTS

  dev@ensembl.org
  mongin@ebi.ac.uk

=cut

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $dbhost      = '';
my $dbuser    = '';
my $dbname    = '';
my $dbpass    = '';
my $dbport;

my $file;
my %band2loc;

&GetOptions (
	    'db:s' => \$dbname,
	    'dbhost:s'=> \$dbhost,
	    'pass:s' => \$dbpass,
            'port:s' => \$dbport,
            'dbuser:s'=> \$dbuser,
	    'file:s'=>\$file,
		    );

print STDERR "Connecting to $dbhost, $dbname\n";


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-port'   => $dbport,
					   );


my %chr_map;
my %map;

open (IN, $file) || die "Can't open $file\n";

my $q1 = "select chromosome_id, name from chromosome";
my $sth = $db->prepare($q1);
$sth->execute();

while (my @ar = $sth->fetchrow_array) {
    $chr_map{$ar[1]} = $ar[0];
}

while (<IN>) {
    chomp;
    my ($chr,$tmpid,$band_start,$band_end,$ori) = split;

    
    my $q2 = "select distinct superctg_name from assembly where chromosome_id != 6 and superctg_name like '%$tmpid'";

    my $sth2 = $db->prepare($q2);
    $sth2->execute();
    
    my $id = $sth2->fetchrow_array();
    
    my $chr_id = $chr_map{$chr};
    if (! defined $chr_id) {
	die;
    }

#Karyotype table: chromosome_id chr_start chr_end band stain
#Get the start and end of the band in chromosome coordinates
    my $q3 = "select min(chr_start),max(chr_end) from assembly where superctg_name = '$id'";
    
    my $sth3 = $db->prepare($q3);
    $sth3->execute();

    my ($min,$max) = $sth3->fetchrow_array();

    my $mid = $min + (($max - $min)/2);

    my $f_start = $min;
    my $f_end = $mid;
    my $s_start = $mid + 1;
    my $s_end  = $max;
    

    if ((! defined $band2loc{$band_start})) {
	$band2loc{$band_start}->{'start'} = $f_start;
	$band2loc{$band_start}->{'end'} = $f_end;
	$band2loc{$band_start}->{'chr_id'} = $chr_id;
    }
    elsif ($band2loc{$band_start}->{'end'} <= $f_end) {
	$band2loc{$band_start}->{'end'} = $f_end;
    }
    elsif ($band2loc{$band_start}->{'end'} <= $f_start) {
	$band2loc{$band_start}->{'end'} = $f_start;
    } 

    if ((! defined $band2loc{$band_end})) {
	$band2loc{$band_end}->{'start'} = $s_start;
	$band2loc{$band_end}->{'end'} = $s_end;
	$band2loc{$band_end}->{'chr_id'} = $chr_id;
    }
    elsif ($band2loc{$band_end}->{'end'} <= $s_end) {
	$band2loc{$band_end}->{'end'} = $s_end;
    }
     elsif ($band2loc{$band_start}->{'end'} <= $s_start) {
	 $band2loc{$band_start}->{'end'} = $s_start;
    } 

}

foreach my $k(keys %band2loc) {
    print "$band2loc{$k}->{'chr_id'}\t$band2loc{$k}->{'start'}\t$band2loc{$k}->{'end'}\t$k\t\\N\n";
    
}



#!/usr/local/bin/perl -w

#Wrapper around Ewan's Birney scanwisepep program

use strict;
use Getopt::Std;

BEGIN {
    my $script_dir = $0;
    $script_dir =~ s/(\S+\/)\S+/$1/;
    unshift (@INC, $script_dir);
    require "mapping_conf.pl";
}

my %conf =  %::mapping_conf;

my $sptr_fa   = $conf{'sptr_fa'};
my $refseq_fa = $conf{'refseq_fa'};

my $organism = $conf{'organism'};
my $target = $conf{'query'};

my $query = $conf{'pmatch_input_fa'};
my $t_thr = $conf{'target_idt'};
my $q_thr = $conf{'query_idt'};
my $opt_o = $conf{'pmatch_out'};

my %target2length;
my %query2length;

$t_thr = 40;
$q_thr = 40;
$target = "/Users/emmanuelmongin/work/test_mapping/Primary/final_test.fa";
$query = "/Users/emmanuelmongin/work/test_mapping/Primary/sptr_ano_gambiae_19_11_02_formated.fa";
my $bin = "/Users/emmanuelmongin/src/wise2.2.3-rc6/src/bin/scanwisep";
$opt_o = "/Users/emmanuelmongin/work/test_mapping/Output/scanpep.out";


#################################
# run scanwisepep

    my $scanwise = "$bin -seqdb $target $query -seqloadtile 5 -hspthread -hspthreadno 4 -hsp2hit_best -hsp2hit_best_perc 10 -hitoutput tab  >> /tmp/$$.scanpep";

print STDERR "Running Scanwise: $scanwise\n";

#system "$scanwise";

#open (PMATCH , "/tmp/$$.scanpep") || die "cannot read /tmp/$$.scanpep\n";

open (OUT,">$opt_o") || die "cannot open $opt_o\n";


open (PMATCH , "/tmp/475.scanpep");

print STDERR "Parsing Output\n";

my $prev_q;
my $prev_t;

my $qlength;
my $tlength;
my $count = 0;

my %match2desc;

while (<PMATCH>) {
    print $_;
    $count++;
    chomp;

    my @a = split;
    my $score = $a[0];
    my $qid = $a[1];
    my $qstart = $a[2];
    my $qend = $a[3];
    my $qlenght = $a[5];
    my $tid = $a[6];
    my $tstart = $a[7];
    my $tend = $a[8];
    my $tlenght = $a[10];
    my $status = $a[11];

    my @matches = ("$qid","$tid");
    
    my $match = join(":",@matches);

    $match2desc{$match}->{'tlength'} = $tlength;
    
    
    $match2desc{$match}->{'qlength'} = $qlength;
    
    my $tmatchlength = $tend - $tstart + 1;
    my $qmatchlength = $qend - $qstart + 1;

    $match2desc{$match}->{'tmatchlength'} += $tmatchlength;
    $match2desc{$match}->{'qmatchlength'} += $qmatchlength;
    $match2desc{$match}->{'status'} = $status;
}

foreach my $key(keys %match2desc) {
    print STDERR "KEY: $key\t LENGHT: $match2desc{$key}->{'tlength'}\n";
    
    my ($query,$target) = split(/:/,$key);
    my $targetperc =   $match2desc{$key}->{'tmatchlength'} * 100 / $match2desc{$key}->{'tlength'}; 
    my $queryperc =  $match2desc{$key}->{'qmatchlength'} * 100 / $match2desc{$key}->{'qlength'}; 
    my $status = $match2desc{$key}->{'status'};
    if (($targetperc >= $t_thr)&&($queryperc >= $q_thr)) {
	print OUT "$query\t$target\t$status\t$queryperc\t$targetperc\n";
    }
}
close(OUT);
print STDERR "The scanwise output $opt_o\nQuery: $query\tThr: $q_thr\nTarget: $target\tThr: $t_thr\n";







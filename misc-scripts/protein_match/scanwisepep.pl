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
my $query = $conf{'query'};

my $target = $conf{'pmatch_input_fa'};
my $t_thr = $conf{'target_idt'};
my $q_thr = $conf{'query_idt'};
my $opt_o = $conf{'pmatch_out'};

my %target2length;
my %query2length;

$target = "/acari/work4/mongin/mouse_5.3/mapping/Primary/target_test.fa";
$query = "/acari/work4/mongin/mouse_5.3/mapping/Primary/t.fa";

#################################
# run scanwisepep

    my $scanwise = "/acari/work4/mongin/tmp/scanwisep -seqdb $target $query -seqloadtile 5 -hspthread -hspthreadno 4 -hsp2hit_best -hsp2hit_best_perc 10 -hitoutput tab  >> /acari/work4/mongin/mouse_5.3/mapping/Output/$$.pmatch";

print STDERR "Running Scanwise: $scanwise\n";

system "$scanwise";

open (PMATCH , "/acari/work4/mongin/mouse_5.3/mapping/Output/$$.pmatch") || die "cannot read $$.pmatch\n";

#open (PMATCH , "/acari/work4/mongin/mouse_5.3/mapping/Primary/scanwise_test.out");

print STDERR "Parsing Output\n";

my $prev_q;
my $prev_t;

my $qlength;
my $tlength;
my $count = 0;

my %match2desc;

my $getseqs = "/usr/local/ensembl/bin/getseqs";

while (<PMATCH>) {
    $count++;
    chomp;
    my @a = split;
    my $score = $a[0];
    my $qid = $a[1];
    my $qstart = $a[2];
    my $qend = $a[3];
    my $tid = $a[5];
    my $tstart = $a[6];
    my $tend = $a[7];
    my $status = $a[9];

    my @matches = ("$qid","$tid");
    
    my $match = join(":",@matches);

    if (! defined $match2desc{$match}->{'tlength'}){
	
	my $cmd = "$getseqs '$tid' $target";
	
	 open(IN,"$cmd 2>/dev/null |") || die "Can't get $tid\n";
	my $seqstr;
	 while(<IN>){
	     chomp;
	     $seqstr .= $_;
	 }
	my $length = length($seqstr);
	#print STDERR "LENGTH: $length\n";
	$match2desc{$match}->{'tlength'} = $length;
    }

    
    if (! defined $match2desc{$match}->{'qlength'}){
	my $cmd = "$getseqs '$qid' $query";

	open(IN,"$cmd 2>/dev/null |") || die "Can't get $qid\n";
	 my $seqstr;
	while(<IN>){
	    chomp;
	    $seqstr .= $_;
	}
	my $length = length($seqstr);
	
	$match2desc{$match}->{'qlength'} = $length;
    }

    my $matchlength = $tend - $tstart + 1;
        
    $match2desc{$match}->{'matchlength'} += $matchlength;
    $match2desc{$match}->{'status'} = $status;
}

foreach my $key(keys %match2desc) {
    print STDERR "$key\t";

    my ($query,$target) = split($key,':');

    my $targetperc =   $match2desc{$key}->{'matchlength'} * 100 / $match2desc{$key}->{'tlength'}; 
    
    my $queryperc =  $match2desc{$key}->{'matchlength'} * 100 / $match2desc{$key}->{'qlength'}; 

    my $st = $match2desc{$key}->{'status'};

    print STDERR "$query\t$queryperc\t$target\t$targetperc\t$status\n";

}








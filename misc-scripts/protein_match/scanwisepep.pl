#!/usr/local/bin/perl -w

#Wrapper around Ewan's Birney scanwisepep program

use strict;
use Getopt::Std;
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;

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

#$t_thr = 40;
#$q_thr = 40;
#$target = "/acari/work4/mongin/final_build/release_mapping/Primary/final.fa";
#$query = "/acari/work4/mongin/final_build/release_mapping/Primary/sptr_ano_gambiae_19_11_02_formated.fa";

#################################
# run scanwisepep

    my $scanwise = "/usr/local/ensembl/bin/scanwisep-2.2.3-rc6 -seqdb $target $query -seqloadtile 5 -hspthread -hspthreadno 4 -hsp2hit_best -hsp2hit_best_perc 10 -hitoutput tab  >> /tmp/$$.pmatch";

print STDERR "Running Scanwise: $scanwise\n";

system "$scanwise";

open (PMATCH , "/tmp/$$.pmatch") || die "cannot read /tmp/$$.pmatch\n";

open (OUT,">$opt_o") || die "cannot open $opt_o\n";


#open (PMATCH , "/tmp/4367756.pmatch");

print STDERR "Parsing Output\n";

my $prev_q;
my $prev_t;

my $qlength;
my $tlength;
my $count = 0;

my %match2desc;

my $format = "fasta";
my @targetdb = "/acari/work4/mongin/final_build/release_mapping/pred_index";
my @querydb = "/acari/work4/mongin/final_build/release_mapping/index";

    my $tfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(
										-db     => \@targetdb,
										-format => $format,
										);

    my $qfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(
										-db     => \@querydb,
										-format => $format,
										);
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

	my $seq = $tfetcher->get_Seq_by_acc($tid);
	my $length = length($seq->seq);
	
	$match2desc{$match}->{'tlength'} = $length;

    }
    
    if (! defined $match2desc{$match}->{'qlength'}){
	my $seq = $qfetcher->get_Seq_by_acc($qid);

	my $length = length($seq->seq);
	
	$match2desc{$match}->{'qlength'} = $length;
    }

    my $tmatchlength = $tend - $tstart + 1;
    my $qmatchlength = $qend - $qstart + 1;

    $match2desc{$match}->{'tmatchlength'} += $tmatchlength;
    $match2desc{$match}->{'qmatchlength'} += $qmatchlength;
    $match2desc{$match}->{'status'} = $status;
}

foreach my $key(keys %match2desc) {
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







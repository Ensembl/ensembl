#!/usr/local/bin/perl

# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself

# wrapper around Richard Durbin's pmatch code (fast protein matcher).

# 06/01 Adapted for Ensembl (mongin@ebi.ac.uk)



use strict;
use Getopt::Long;
#use vars qw($opt_q $opt_t $opt_l $opt_o $t_thr $q_thr $opt_w $opt_s $opt_c $opt_d);

BEGIN {
    my $script_dir = $0;
    $script_dir =~ s/(\S+\/)\S+/$1/;
    unshift (@INC, $script_dir);
    require "mapping_conf.pl";
}



my %conf =  %::mapping_conf;

my $sptr_fa   = $conf{'sptr_fa'};
my $refseq_fa = $conf{'refseq_fa'};


#Set the default percentage of idt

my $opt_q = $conf{'ensembl_predictions'};
my $opt_t = $conf{'total_known_fa'};
my $opt_o = $conf{'mapping_out'};

my $t_thr = $conf{'min_ensembl_idt'};
my $q_thr = $conf{'min_known_idt'};
my $pmatch_bin = $conf{'pmatch'};
my ($opt_w,$opt_l,$opt_d);
my $help;


&GetOptions(
	    'help' => \$help,
	    );

if ($help) {
    print STDERR $conf{'help'}."\n";
    exit();
}

#Check if the configuration file is correct
my %check;
$check{'ensembl_predictions'} = $conf{'ensembl_predictions'};
$check{'total_known_fa'} = $conf{'total_known_fa'};
$check{'mapping_out'} = $conf{'mapping_out'};
$check{'target_id'} = $conf{'min_ensembl_idt'};
$check{'min_known_idt'} = $conf{'min_known_idt'};
$check{'pmatch'} = $conf{'pmatch'};

foreach my $k (keys %check) {
    if ($check{$k} !~ /(\S+)/) {
	usage();
    }
}


#End of checks

#################################

my $query = $opt_q;
my $target = $opt_t;
my %hash2;


#################################
# run pmatch (Richard Durbin's fast protein matcher, rd@sanger.ac.uk)

print STDERR "run pmatch...\n";
if ($opt_l) {
    my @files;
    open (DB , "$opt_l") || die "cannot read $opt_l\n";
    while (<DB>) {
        chomp;
        my @a = split /\s+/;
        push (@files , @a);
    }
    close DB;        
    foreach my $file (@files) {
        my $pmatch = $pmatch_bin."  -T 14 $file $query >> $$.pmatch";
	print $pmatch."\n";
        system "$pmatch";
    }
}
else {
    my $pmatch = $pmatch_bin."  -T 14 $target $query > $$.pmatch";
    print $pmatch."\n";
    system "$pmatch";
}

if ($opt_w && !$opt_d) {
    unlink "$$.swall";
}

#################################
# sort the pmatch output file, based on query and target name
# (necessary since not all the matches relating to a query-target pair cluster)

print STDERR "sort the pmatch output file...\n";
my (@a , @q , @t);
open (TMP , ">$$.sort") || die "cannot create $$.sort\n";
open (PMATCH , "$$.pmatch") || die "cannot read $$.pmatch\n";
while (<PMATCH>) {
    my @f = split /\t/;
    push @a, $_;
    push @q, $f[1];
    push @t, $f[6];
}
foreach my $i (sort { $q[$a] cmp $q[$b] or $t[$a] cmp $t[$b] } 0..$#a) { print TMP $a[$i] }
close PMATCH;
close TMP;
rename ("$$.sort" , "$$.pmatch");

#################################
# determine the best path of non-overlapping matching substrings 

print STDERR "stitch pmatch substrings...\n";
open (STITCH , ">$$.stitch") || die "cannot create $$.stitch\n";
my @match_list = ();
my $old_query;
my $old_target;
my $old_qlen;
my $old_tlen;
my $previous_line;
open (PMATCH , "$$.pmatch") || die "cannot read $$.pmatch\n";
while (<PMATCH>) {
    chomp;
    my @a = split /\t/;
    my $query = $a[1];
    my $qstart = $a[2];
    my $qend = $a[3];
    my $qlen = $a[5];
    my $target = $a[6];
    my $tstart = $a[7];
    my $tend = $a[8];
    my $tlen = $a[10];
    
    # new set of query/target 
    if (($query ne $old_query || $target ne $old_target) && @match_list) {
	if (@match_list == 1) {
	    print STITCH "$previous_line\n";
	}
	else {
	    my ($max , $trace) = stitch_matches (@match_list);
	    my $qperc = sprintf ("%.1f" , ($max/$old_qlen)*100);
	    my $tperc = sprintf ("%.1f" , ($max/$old_tlen)*100);
	    my $num = @$trace-1;
	    print STITCH "$max\t$old_query\t$trace->[1]->{QSTART}\t$trace->[$num]->{QEND}\t$qperc\t$old_qlen\t";
	    print STITCH "$old_target\t$trace->[1]->{TSTART}\t$trace->[$num]->{TEND}\t$tperc\t$old_tlen\t($num)\n";
	    }
	@match_list = ();
	my $match = Pmatch->new('QSTART'=>$qstart,'QEND'=>$qend,'TSTART'=>$tstart,'TEND'=>$tend);
	push (@match_list , $match);
	$old_query = $query;
	$old_target = $target;
	$old_qlen = $qlen;
	$old_tlen = $tlen;
	$previous_line = $_;
    }
    # last line 
    elsif (eof) {
	unless (@match_list) {
	    $old_query = $query;
	    $old_target = $target;
	    $old_qlen = $qlen;
	    $old_tlen = $tlen;
	}
	my $match = Pmatch->new('QSTART'=>$qstart,'QEND'=>$qend,'TSTART'=>$tstart,'TEND'=>$tend);
	push (@match_list , $match);
	my ($max , $trace) = stitch_matches (@match_list);
	my $num = @$trace-1;
	my $qperc = sprintf ("%.1f" , ($max/$old_qlen)*100);
	my $tperc = sprintf ("%.1f" , ($max/$old_tlen)*100);
	print STITCH "$max\t$old_query\t$trace->[1]->{QSTART}\t$trace->[$num]->{QEND}\t$qperc\t$old_qlen\t"; 
	print STITCH "$old_target\t$trace->[1]->{TSTART}\t$trace->[$num]->{TEND}\t$tperc\t$old_tlen\t($num)\n";
    }
        # else
    else {
	my $match = Pmatch->new('QSTART'=>$qstart,'QEND'=>$qend,'TSTART'=>$tstart,'TEND'=>$tend);
	push (@match_list , $match);
	$old_query = $query;
            $old_target = $target;
	$old_qlen = $qlen;
	$old_tlen = $tlen;
	$previous_line = $_;
    }
}
close PMATCH;
close STITCH;


##########################################
# process the matches, classifying them into several categories


print STDERR "process matches...\n";
# get first a list of all query id's
my @queries;
    open (ID , "$query") || die "cannot read $query\n";
while (<ID>) {
    chomp;
    if (/^\>(\S+)/) {
	push (@queries , $1);
    }
}
close ID;

# process the matches
my $input;
#if ($opt_s) {
$input = "$$.stitch";
    #}
#else {
#   $input = "$$.pmatch";
#}
open (PROCESS , ">$$.process") || die "cannot create $$.process\n";
open (IN , "$input") || die "cannot read $input\n";
process_matches (*IN , *PROCESS , \@queries);
close IN;
close PROCESS;

#########################################
# delete some tmp files, and move the appropriate results to the out file unless ($opt_d)

rename ("$$.process" , "$opt_o");
unlink "$$.pmatch";
unlink "$$.stitch";

rename ("$$.stitch" , "$opt_o");
unlink "$$.pmatch";        

########################
# subroutines
########################
sub stitch_matches {
    my @match_list = @_;
    my $DEBUG = 0;
    my $TRACE = 1;

    # sort the matches, based on the start coordinate, and keep them in the @sort array
    # !! does not explicitly deal with cases where match boundaries overlap
    #    (like qstart1 < qstart2 and tstart1 > tstart2). But since we only
    #    chain non-overlapping matches, this should be ok (the order of the matches
    #    within the clusters of overlapping matches does not matter)
    my @sort = ();
    foreach my $match (@match_list) {
        unless (@sort) {
            push (@sort , $match_list[0]);
            next;
        }
        my $element_num = @sort;
        my $switch = 0;
        for (my $i = 0 ; $i < $element_num ; $i++) {
            my $sorted_match = $sort[$i];
            if ($match->{QSTART} < $sorted_match->{QSTART} && $match->{TSTART} < $sorted_match->{TSTART}) {
                $switch = 1;
                # insert the new match into the list (kind of linked list)
                my @tmp = splice (@sort , $i);
                push (@sort , $match);
                push (@sort , @tmp);
                last;
            }
        }
        unless ($switch) {
            push (@sort , $match);
            $switch = 0;
	}
    }
    
    if ($DEBUG) {
        foreach (@sort) {
            print STITCH "\tsort: $_->{QSTART}  $_->{QEND}  $_->{TSTART}  $_->{TEND}\n";
	}
    }

    # loop over all the matches, always calculating the best score up to this point
    # (some kind of dynamic programming):
    #     - accept only non-overlaping matches
    #     - use the length of the uniquely matching sequence as score
    #       (not really necessary with pmatch, since it only returns exact matches)
    #       # my $q = $sort[$j]->{QEND} - $sort[$j]->{QSTART} + 1;
    #       # my $t = $sort[$j]->{TEND} - $sort[$j]->{TSTART} + 1;
    #       # my $tmp_score = $score[$i] + (($q < $t) ? $q : $t);
    #       my $tmp_score = $score[$i] + $sort[$j]->{QEND} - $sort[$j]->{QSTART} + 1;
    my $max = 0;
    my $max_index = 0;
    my @score;
    my @trace;
    my $tmp_trace;

    # define the (starting) boundary condition, and make this match the first element of @sort
    my $boundary = Pmatch->new('QSTART'=>0,'QEND'=>0,'TSTART'=>0,'TEND'=>0);
    unshift (@sort , $boundary);
    $score[0] = 0;
    push (@{$trace[0]} , $sort[0]);

    # loop over all the matches, always calculating the best score up to this point
    # (some kind of dynamic programming)
    for (my $j = 1 ; $j < @sort ; $j++) {
        $score[$j] = 0;
        for (my $i = 0 ; $i < $j ; $i++) {
            # accept only non-overlaping matches
            if ($sort[$i]->{QEND} < $sort[$j]->{QSTART} && $sort[$i]->{TEND} < $sort[$j]->{TSTART}) {
                my $tmp_score = $score[$i] + $sort[$j]->{QEND} - $sort[$j]->{QSTART} + 1;
                # keep the best score, and the trace pointer
                if ($tmp_score > $score[$j]) {
                    $score[$j] = $tmp_score;
                    $tmp_trace = $i;
                }
            }
        }
        # make the trace update
        if ($TRACE) {
            push (@{$trace[$j]} , @{$trace[$tmp_trace]});
            push (@{$trace[$j]} , $sort[$j]);
	}
        # recalculate the max score
        if ($score[$j] > $max) {
            $max = $score[$j];
            $max_index = $j;
        }
    }

    if ($DEBUG) {
        for (my $j = 1 ; $j < @score ; $j++) {
            print STITCH "match $j (score $score[$j]):\n";
            for (my $i = 1 ; $i < @{$trace[$j]} ; $i++) {
                print STITCH "\ttrace $trace[$j]->[$i]->{QSTART}  $trace[$j]->[$i]->{QEND}  $trace[$j]->[$i]->{TSTART}  $trace[$j]->[$i]->{TEND}\n";
	    }
        }
    }
    if ($TRACE) {
        return ($max , $trace[$max_index]);    
    }
    else {
        return ($max);
    }
}

####################
sub process_matches {
    local *IN = shift;
    local *OUT = shift;
    my $query_ids = shift;

    # define some variables:
    # to keep track of the quality of the matches
    my %match;
    my %partial;
    my %candidate;
    my %target_match;

    # for the parsing
    my %percent;
    my $old_qid = "";
    my $old_tid = "";
    my $best_qperc = 0; 

    # to count how often the different match classes occurred
    my $count_match = 0;
    my $count_partial_match = 0;
    my $count_partial = 0;
    my $count_candidate = 0;
    my $count_orphan = 0;

    # read the pmatch output file,
    # and keep the "percent matching" in a 2D hash, indexed by query and target sequence
    while (<IN>) {
        chomp;
        my ($len,$qid,$qstart,$qend,$qperc,$qlen,$tid,$tstart,$tend,$tperc,$tlen) = split /\t/;
        # we only consider the best match per query-target pair. Use the -s option to
        # get the total matching sequence (pmatchn gives multiple hits either due to e.g.
        # introns, or to internal repeats (we don't want the latter))
        if ($qid ne $old_qid || $tid ne $old_tid) {
            $percent{$qid}->{$tid} = $qperc."\t".$tperc;
            $best_qperc = $qperc;
	}
        elsif ($qid eq $old_qid && $tid eq $old_tid && $qperc > $best_qperc) {
            $percent{$qid}->{$tid} = $qperc."\t".$tperc;
            $best_qperc = $qperc;
	}
        $old_qid = $qid;
        $old_tid = $tid;
    }

    # classify the matches
    foreach my $query_id (sort {$a cmp $b} keys %percent) {
        my @matches = ();
        my @partial_query = ();
        my @partial_target = ();
        my @candidates = ();
        my @best = ();   
        # loop over all targets, having the matches sorted (descending) based on the query percent coverage ($qperc)
        # and populate the different arrays of types of matches
        foreach my $target_id (sort {$percent{$query_id}{$b} <=> $percent{$query_id}{$a}} keys %{$percent{$query_id}}) {
            my ($qperc , $tperc) = split (/\t/ , $percent{$query_id}->{$target_id});

            if ($qperc > $q_thr && $tperc > $t_thr) {
                
		if( !defined $hash2{$target_id} ) {
		    $hash2{$target_id} = [];
		}
		
		my $p= NamePerc->new;
		$p->name($query_id);
		$p->tperc($tperc);
		$p->qperc($qperc);

		
		push(@{$hash2{$target_id}},$p);

	    }
            else {
                last;
	    }
	}
    }

#Loop over the target protein
    foreach my $know ( keys %hash2 ) {
	my @array = @{$hash2{$know}};
	@array = sort { $b->tperc <=> $a->tperc } @array;
	
#The Ensembl match to the known protein is labelled as PRIMARY and will be used later for the mapping 
	my $top = shift @array;
	
#	print  OUT $top->name."\t$know\tcandidate\t".$top->qperc."\t".$top->tperc."\n";
	print OUT "$know\t".$top->name."\tcandidate\t".$top->tperc."\t".$top->qperc."\n";

	foreach my $repeat (@array) {
	    if( $repeat->tperc+1 >= $top->tperc ) {
#		print  OUT $repeat->name."\t$know\tduplicate\t".$repeat->qperc."\t".$repeat->tperc."\n";
		print OUT "$know\t".$repeat->name."\tduplicate\t".$repeat->tperc."\t".$repeat->qperc."\n";
	    }
	}
    }
}

sub usage {
    
  print STDERR <<HELP

Usage: pmatch.pl 
One of the element of the configuration file has not been properly loaded
Please fill in properly your configuration file

Here is your set up:
HELP
;

 foreach my $k (keys %check) {
	print STDERR "$k:\t$check{$k}\n";
    }



  exit();
}


##########################################
# some OO stuff:
##########################################

package Pmatch;

# constructor
sub new {
    my $class = shift;
    my $self = {};
    bless ($self , $class);
    $self->_init(@_);  # call _init to initialise some attributes
                       # interpret the remaining args as key-value pairs
    return $self;
}

# _init method, initialising some attributes,
# and interpreting the remaining args as key-value pairs

sub _init {
    my $self = shift;
    my %extra = @_;
    @$self{keys %extra} = values %extra;
}


package NamePerc;


sub new {
    my $class= shift;
    my $self = {};
    bless $self,$class;
    return $self;
}


=head2 name

 Title   : name
 Usage   : $obj->name($newval)
 Function:
 Returns : value of name
 Args    : newvalue (optional)


=cut

sub name{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'name'} = $value;
    }
    return $obj->{'name'};

}

=head2 qperc

 Title   : qperc
 Usage   : $obj->qperc($newval)
 Function: 
 Returns : value of query perc
 Args    : newvalue (optional)


=cut

sub qperc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'qperc'} = $value;
    }
    return $obj->{'qperc'};

}


=head2 tperc

 Title   : tperc
 Usage   : $obj->tperc($newval)
 Function:
 Returns : value of target perc
 Args    : newvalue (optional)


=cut

sub tperc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'perc'} = $value;
    }
    return $obj->{'perc'};

}                    


1;
##########################################

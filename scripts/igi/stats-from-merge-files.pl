#!/usr/local/bin/perl
# $Id$ 

# Copyright EnsEMBL http://www.ensembl.org

# This script is used for collecting statistics on the gtf merges done by
# gtf_merge.pl (which are typically called from the wrapper
# all-merges.sh), and it is typically called by the all-stats.sh wrapper.

# The whole objective is to compare gene predictions of different
# 'sources', currently EnsEMBL, Affymetrix and Fgenesh (Softberry Inc),
# and arrive at a consensus 'Initial Gene Index' (igi). The igi's are
# assigned by gtf_merge, and this script just summarizes them. 

# The script has grown a bit unwieldy, but that is because it is so
# convenient to stay inside this script once every feature, gene and igi
# of each prediction is loaded into the hashes-of-hashes-of-hashes of this
# script.

# For usage, ask author ;-)

# Written by Philip lijnzaad@ebi.ac.uk


use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::igi_utils;


my $usage = "$0 options < merged-file.gtf  > outputfile\n";
my $help;
my $all_stats;  
my $chaining = undef;
my $i2nmapping;
my $n2imapping;
my $gtf_summary;
my $cluster_n;
&GetOptions( 
            'stats'  => \$all_stats,
            'chaining:s'  => \$chaining,
             'gtfsummary:s' => \$gtf_summary,
             'cluster_n:s'  => \$cluster_n,
	     'igi2native:s'     => \$i2nmapping, 
	     'native2igi:s'     => \$n2imapping, 
	     'h|help'     => \$help,
	     );
die $usage if $help;

die "need both cluster_n and gtf_summary or neither"
  unless ( defined($gtf_summary) == defined($cluster_n));

die "need -cluster_n when doing igi2native or native2igi"
  if ( (defined($n2imapping)||defined($i2nmapping)) >   
           defined($cluster_n));


my @argv_copy = @ARGV; # may get gobbled up by the <> construct. 

my %igis_of_source; # $igis_of_source_{$source}{$igi} = [nfeats, start,
                    # end, nexons ]
my %natives_of_source; # $natives_of_source{$source}{$native_id} =
                       # [ nfeats, start, end, nexons ]
my %natives_of_igi; # $natives_of_igi{$source}{$igi}{$native_id} => 1

my %transcripts_of_native; #
                           # $transcripts_of_native{$source}{$native}{$transcript}
                           # = length of the transcript

my $blurp = blurp();

### just read the complete file into the relevant hashes
GTF_LINE:
while (<>) {
    # taken largely from Bio::EnsEMBL::Utils::GTF_handler:
    next GTF_LINE if /^#/;
    next GTF_LINE if /^\s*$/;
    chomp;
    
    my @fields = split "\t", $_;
    my ($seq_name, $source, $feature,
        $start,  $end,    $score,
        $strand, $phase,  $group_field)  = @fields;
    $feature = lc $feature;

    unless ($group_field) {
        warn("no group field: skipping : '$_'\n");
        next GTF_LINE ;
    }
    
    # Extract the extra information from the final field of the GTF line.
    my ($igi, $gene_name, $native_ids, $transcript_ids, $exon_num, $exon_id) =
      Bio::EnsEMBL::Utils::igi_utils::parse_group_field($group_field);

    unless ($igi) {
        warn("Skipping line with no igi_id: '$_'\n");
        next GTF_LINE;
    }

    my $native_id; 
    if (@$native_ids == 0 ) {
        warn("Skipping line with no gene_id(s): '$_'\n");
        next GTF_LINE;
    }

    if (int(@$native_ids) > 1 ) {
        warn("Line with several gene_ids (taking first one): '$_'\n");
    }
    $native_id =  $ {$native_ids}[0];

    if  (@$transcript_ids == 0) {
        warn("Skipping line with no transcript_id(s): '$_'\n");
        next GTF_LINE;
    }
    if (int(@$transcript_ids) > 1 ) {
        warn("Line with several transcript_ids (taking first one): '$_'\n");
    }
    my $transcript_id =  $ {$transcript_ids}[0];

    # keep track of marginals by 'looping' over current one +
    # fictional source 'ALL':
    foreach my $s ($source, 'ALL') {

        # keep track of number of features, exons, start and end of this igi:
        track_extents(\%igis_of_source, $s, $igi,
                      $seq_name, $start, $end, $strand, $exon_num);

        # as well as this native_id:
        track_extents(\%natives_of_source, $s, $native_id,
                      $seq_name, $start, $end, $strand, $exon_num);

        # keep track of the native gene_ids of a given igi in a given
        # source. As always, use a hash for faster collating (i.e. keys
        # %($natives_of_igi{$s}{$igi}) gives them all; the value is
        # irrelevant 
        $natives_of_igi{$s}{$igi}{$native_id}++;
        
        # pointless to keep track of exon statistics; call exon-lengths.awk
        # for that.
    }                                   # foreach $s

    if ( $feature eq 'exon' ) {
        $transcripts_of_native{$source}{$native_id}{$transcript_id} +=
          ($end - $start);
    }
}                                       # while(<>)
### done reading files

# get rid of 'ALL' (which was added for convenience when gathering stats
# on marginals)
my @all_sources = (grep $_ ne 'ALL', sort keys %igis_of_source);

my %sources_of_igi= undef;
my @igis_of_n_sources = undef;
&find_overlaps;

my @native_ids_per_igi_histo = undef;
my @igis_per_chaining_group = undef;
&find_chaining;

if ($all_stats) {
    print_stats();
}

if ($chaining) {
    print "igi's and native ids of igi's that chain together >= $chaining native ids\n";
    # print out native_id's clumped together in clusters of more than n
    foreach my $source ( @all_sources ) {
        print_chaining($source, $chaining);
    }
    print "----\n";
}

if ($i2nmapping) {
    foreach my $s (@all_sources) { 
        my $file = "> $i2nmapping.$s";
        open(I2N, $file) || die "can't open $file: $!";
        print_igi_to_native(\*I2N, $s, $cluster_n);
        close(I2N);
    }
}

if ($n2imapping) {
    foreach my $s (@all_sources) {
        my $file = "> $n2imapping.$s";
        open(N2I, $file) || die "can't open $file: $!";
        print_native_to_igi(\*N2I, $s, $cluster_n);
        close(N2I);
    }
}

if ($gtf_summary) {
    my $file = "$gtf_summary";
    open(OUT, "> $file" ) || die die "$file: $!";
    
    my (@stuff)  = ("on", `date`, "by", `whoami`, "@",  `hostname`);
    foreach (@stuff) {chomp};
    print OUT $blurp;
    print OUT <<MOREBLURP
## 
## GTF summary file in gff format. The source is 'igi3', and only gene
## features are given. Only genes predicted by $cluster_n or more gene
## predictions have been included. The original gene id's have been
## prefixed with <source-name>. The output is sorted by fpc contig id,
## strand, and start.  This is not the final IGI3 file; it just summarizes
## the extent of the igi3 clusters.
MOREBLURP
  ;
    close(OUT);
    my $sortcmd="sort -k1,1 -k7,7 -k4,4n ";
    open(OUT, "| $sortcmd >> $file") || die "$file: $!";
    gtf_for_igis_predicted_by_n(\*OUT, $cluster_n);
}                                       # if gtf_summary

### simple log message
sub blurp {
    my (@stuff)  = ("on", `date`, "by", `whoami`, "@",  `hostname`);
    foreach (@stuff) {chomp};
    my $s =  '### Produced by $Id$  ' . "\n" 
      . "### run " . join(' ',@stuff) .  "\n"
      . "### for EnsEMBL (http://www.ensembl.org)\n"
      . "### argument(s): ". join(' ', @argv_copy) . "\n";
    foreach (@argv_copy) {   $s .= "### ", `ls -l $_`; }
    $s .= "###\n";
    $s;
}

sub print_coord_stats {
    my ($min, $max, $avg, 
        $minfeats, $maxfeats, $avgfeats, 
        $minexons, $maxexons, $avgexons)  = @_;

    print "\tminl=$min, maxl=$max, avgl=$avg\n";
    print "\tminfeats=$minfeats, maxfeats=$maxfeats, avgfeats=$avgfeats\n";
    print "\tminexons=$minexons, maxexons=$maxexons, avgexons=$avgexons\n";
}

sub igi_stats {
    my ($igi_hash) = @_;
    my $huge = 1e308; 

    my ( @igis) = keys %$igi_hash;
    my $min = $huge;
    my $max = -$huge;
    my $sum = 0;
    my $minfeats = $huge;
    my $maxfeats = -$huge;
    my $sumfeats = 0;

    my $minexons = $huge;
    my $maxexons = -$huge;
    my $sumexons = 0;

    my $n =0; 
    foreach my $igi (@igis) {
        my ($nfeats, $start, $end, $nexons) = @{$igi_hash->{$igi}};
        die unless defined($nexons);
        my $len = ($end - $start);
        $min = $len if $len <  $min ;
        $max = $len if $len >  $max ;
        $sum += $len;

        $minfeats = $nfeats if $nfeats <  $minfeats ;
        $maxfeats = $nfeats if $nfeats >  $maxfeats ;
        $sumfeats += $nfeats;

        $minexons = $nexons if $nexons <  $minexons ;
        $maxexons = $nexons if $nexons >  $maxexons ;
        $sumexons += $nexons;

        $n++;
    }

    my ($avg, $avgfeats, $avgexons) = ('none', 'none', 'none');
    if ($n>0) {
        $avg = int($sum/$n);
        $avgfeats= int($sumfeats/$n);
        $avgexons = int($sumexons/$n)
    } 

    print_coord_stats( $min, $max, $avg, 
                 $minfeats, $maxfeats, $avgfeats, 
                 $minexons, $maxexons, $avgexons);
}


### keep track of start,end of a gene, by looking at the lowest start and
### higest end of any of the features. This is used both for igi's and
### native id's
sub track_extents {
    my ($extents, $source, $id, 
        $seq_name, $start, $end, $strand, $exon_num) = @_;
    my ($nfeats, $min, $max, $nexons);
    
    #get previous record of this igi, if any:
    if (defined $ {$extents}{$source}{$id}) {
        ($nfeats, $min, $max, $nexons) = 
          @{$ {$extents}{$source}{$id}}[0,1,2,3];
    } else { 
        ($nfeats, $min, $max, $nexons)= (0, $start, $end, 0);
    }
        
    $min = $start if $start < $min;
    $max = $end if $end > $max;
    $nfeats++;
    $nexons = $exon_num if $exon_num > $nexons;
    
    # add record back in:
    $ {$extents}{$source}{$id} = [$nfeats, $min, $max, $nexons, 
                                 $seq_name, $strand ];
}

### find out how the sets overlap between each other
sub find_overlaps { 
    my @all_igis = keys %{$igis_of_source{'ALL'}};
    my $n_igis = int(@all_igis);

    # invert igis_of_source:
    foreach my $igi (@all_igis) {
        foreach my $source (@all_sources ) { # excluding 'ALL', of course
            # does this source contain this igi?
            if (defined $ {$igis_of_source{$source}}{$igi} ) {
                $sources_of_igi{$igi}{$source}++; # meaning: seen it
            }
        }
    }

## do histogram: find the igis per overlap group (i.e., n sources)
    foreach my $igi (@all_igis) {
        my $n  = int(keys %{$sources_of_igi{$igi}});
        my $loc = $ {$igis_of_source{'ALL'}}{$igi}; # start,end etc. of this igi
        
        #  $ {$igis_of_n_sources[$n]}{$igi} = $loc 
        # wrong, this is for those that are in
        # exactly two groups, but you want to know the  cumulation:
        for (my $i=$n; $i>=0; $i--) {
            $ {$igis_of_n_sources[$i]}{$igi} = $loc
        }
    }
}                                       # find_overlaps

sub find_chaining { 
    
    foreach my $source ('ALL', @all_sources) { 
        foreach my $igi ( keys % {$igis_of_source{$source}}  ) {
            my @native_ids = keys %{$natives_of_igi{$source}{$igi}};
            # these are all distinct native id's of this igi in this source;
            # collate them in a histogram
            my $n = int(@native_ids);
            $native_ids_per_igi_histo[ $n ]{$source} ++;
            
            my $loc = $ {$igis_of_source{'ALL'}}{$igi}; # start,end etc. of this igi
            $ {igis_per_chaining_group[ $n ]{$source}}{$igi} = $loc;
        }
    }
}

sub print_stats {
    print $blurp;                       # log comment
    print "Sources: " , join( ' ', @all_sources), "\n";

    my $n_sources = int(@all_sources);

    print "number of igi's per source:\n";
    foreach my $source ('ALL', @all_sources) {
        my $igi_of_source = $igis_of_source{$source};
        my $num = int(keys %{$igi_of_source});
        print "$source: $num\n";
    }
    print "----\n";

    # compare igi's per source
    print "number of igis shared with other source; sharing multiplicity\n";
    my %h = undef;
    # headers:
    print '  '
      # overlaps with named groups
      , map(do {sprintf "%12s  ", substr($_,0,6) }, (@all_sources))
        # multiplicity of groups membership:
        , map(do {sprintf "        >=%2d  ", 2+$_ }
              ,    reverse ( 0 .. (int(@all_sources)-2)))
          , "    unshared  ",
          , "       total\n";
    
SOURCE1:
    foreach my $source1 ('ALL', @all_sources) {
        my @igis_of_source1 = keys %{$igis_of_source{$source1}};
        my $source2;
        
        printf "%6s  ", $source1;
        SOURCE2:
        my ($overlap, $tot);
        my %is_shared=undef;
        foreach $source2 (@all_sources) {
            ($overlap, $tot) = (0,0);
            foreach my $igi (@igis_of_source1) {
                if (defined $igis_of_source{$source2}{$igi}) {
                    $overlap++;
                    # keep track of how many sources have this particular igi:
                    $is_shared{$igi}++; #  if $source1 ne $source2;
                }
                $tot++;
            }
            if ($source1 eq $source2) {
                print "    --        ";
            } else { 
                printf "%6d (%2d%%)  ", $overlap, 100.0*$overlap/$tot;
            }
        }
        # print out the sharing multiplicity
        my ($from, $to) = (int(@all_sources), 2);
        for(my $i =$from; $i>= $to; $i--) {
            my $num = int(grep( $_ >= $i, values %is_shared));
            printf "%6d (%2d%%)  ", $num, 100.0*$num/$tot;
        }
        my $unshared =
          int(grep( $_ == 1, values %is_shared)); # those seen exactly once
        # are in just one source,
        # so are not shared
        printf "%6d (%2d%%)  %6d\n", $unshared, 100.0*$unshared/$tot, $tot; 
    }                                   # for source1

    print "----\n";

    ### igi statistics per source
    foreach my $source ('ALL', @all_sources) {
        print "igi stats on $source:\n";
        igi_stats($igis_of_source{$source});
    }
    print "----\n";

    ### see to what extent the clusters grow as they are
    ### shared by an increasing number of sources:
    print "igi stats per cluster group:\n";
    for(my $i = 1; $i<=$n_sources; $i++) {
        print "those in $i sources:\n";
        igi_stats($igis_of_n_sources[$i]) ;
    }
    print "----\n";

    ### output the chaining: the extent to which native gene_ids get
    ### chained together by the igi clustering
    print "gene id chaining: histogram of native gene_ids per igi, and igi statistics on this:\n";
    foreach my $source ('ALL', @all_sources) { 
        print "source $source:\n";
        for (my $i=0;  $i< int(@native_ids_per_igi_histo); $i++ ) {
            my $n = $native_ids_per_igi_histo[ $i ]{$source};
            if ($n > 0) { 
                print "\t$i: $n\n";
                igi_stats( $ {igis_per_chaining_group[ $i ]}{$source}) ;
            }
        }
    }
    print "----\n";
}                                       # all_stats

# print out native_id's clumped together in clusters of more than n
sub print_chaining {
    my ($source, $n)=@_;
    
    print "source $source:\n";
    for (my $i=$n; $i< int(@native_ids_per_igi_histo); $i++) {
        my $num = $native_ids_per_igi_histo[ $i ]{$source};
        if ( $num ) { 
            print "    chaining together $i native ids ( $num cases) :\n";
            my $igi_hash = $ {igis_per_chaining_group[ $i ]}{$source};
            foreach my $igi ( keys %$igi_hash) {
                my ($start,$end, $seq, $strand) 
                  = @{$igis_of_source{$source}{$igi}}[1,2,4,5];
                print "      $igi [ $seq\t$start\t$end\t$strand ]:\n";  
                print_natives($source, $igi);
                print "        linked by:\n";
                foreach my $s ( @all_sources ) {
                    next if $s eq $source;
                    print_natives($s, $igi, "  ");
                }
            }
        }
    }
}                                       # print_chaining

# print out the native id's of this igi, with coordinates, and sorted by
# coordinate
sub print_natives {
    my ($source, $igi, $indent) = @_;

    my (@native_coords);
    foreach my $nat ( keys %{$natives_of_igi{$source}{$igi}} ) {
        my ($start,$end) = @{$natives_of_source{$source}{$nat}}[1,2];
        push @native_coords, [$nat, $start, $end];
    }
    my $n; $n++;
    foreach my $nat (sort sort_native @native_coords) { 
        my ($id, $start, $end) = @{$nat};
        print "        $indent$source $id [ $start\t$end ]\n";
    }
}                                       # print natives

# for sorting id's by start/stop coords. This is a bit hairy, to get the
# start exon's first, and end exons last ... 
sub sort_native {
    my ($starta, $enda) = @{$a}[1,2];
    my ($startb, $endb) = @{$b}[1,2];

    my $n;
    $n = $starta <=> $startb;
    return $n if $n;
    $n = $enda <=> $endb;
    return $n if $n;
    return $$a[0] cmp $$b[0];           # i.e., alphabetically
}


## print file that maps igi to native
sub print_igi_to_native {
    my($OUT, $source, $cluster_n) = @_;
    foreach my $igi (sort 
                     { Bio::EnsEMBL::Utils::igi_utils::by_igi_number($a,$b)}
                     keys %{$igis_of_source{$source}}) {

        if (defined $ {$igis_of_n_sources [ $cluster_n ]}{$igi} ) {
            print $OUT "$igi "
              , join(' ', (sort keys %{$natives_of_igi{$source}{$igi}}))
                , "\n";
        } else { warn "i2n: $igi: not in $cluster_n sources\n" }
    }
}                                       # print_igi_to_native

## print file that maps native to igi:
sub print_native_to_igi {
    my($OUT, $source, $cluster_n) = @_;

    my %igis_of_native = undef;         # have to invert first:

    foreach my $igi (sort 
                     { Bio::EnsEMBL::Utils::igi_utils::by_igi_number($a,$b)}
                     keys %{$igis_of_source{$source}}) {
        
        foreach my $nat (keys %{$natives_of_igi{$source}{$igi}}) {
            $igis_of_native{$nat}=$igi;
        }
    }

    foreach my $nat (sort keys %igis_of_native) {
        my $igi = $igis_of_native{$nat};
        if (defined $ {$igis_of_n_sources [ $cluster_n ]}{$igi} ) {        
            print $OUT "$nat $igi\n";
        } else { warn "n2i: $igi: not in $cluster_n sources\n" }
    }
}                                       # print_native_to_igi

### print a gtf file for all igi's that are predicted by N sources. For
### now, just print start and end, nothing else. I.e., don't give any
### native exons or so:
sub gtf_for_igis_predicted_by_n {
    my ($OUT, $n)  = @_;
    my $newsource = 'igi3';
    my $feature='gene';
    my $score = 0;
    my $phase = '.';

    foreach my $igi  (sort 
                      { Bio::EnsEMBL::Utils::igi_utils::by_igi_number($a, $b)}
                      keys %{$igis_of_n_sources[$n]}) {
        my ($nfeats, $min, $max, $nexons, $seq_name, $strand )
          = @{$igis_of_source{'ALL'}{$igi}};
        my @fields =($seq_name, $newsource, $feature, $min,$max, $score, $strand, $phase);

        # add the ids:
        my $rest ="igi_id \"$igi\"; ";
      SOURCE:
        foreach my $source (sort @all_sources) {
            foreach my $nat_id (sort keys %{$natives_of_igi{$source}{$igi}}) {
                $rest .= "gene_id \"$source:$nat_id\"; ";
                
                ## find id of the longest transcript:
                my $max = -1; my $max_tid=undef;
                my @tids = keys %{$transcripts_of_native{$source}{$nat_id}};
# warn "tids:  @tids" if @tids > 2;
                foreach my $tid ( @tids ) { 
                    my $len = $transcripts_of_native{$source}{$nat_id}{$tid};
                    if ($len > $max) {
                        $max_tid = $tid;
                        $max = $len;
                    }
                }
                $rest .= "transcript_id \"$source:$max_tid\"; ";
            }                           # foreach $nat
        }                               # source
        push(@fields, $rest);
        print $OUT join("\t", @fields), "\n"; 
    }                                   # foreach $igi
}                                       # gtf_for_igis_predicted_by_n

#!/usr/local/bin/perl

# $Id$ 
use strict;
use Getopt::Long;

my $usage = "$0 options < merged-file.gtf  > outputfile\n";
my $help;
my $all_stats;  
my $chaining = undef;

&GetOptions( 
            'stats'  => \$all_stats,
            'chaining:s'  => \$chaining,
	     'h|help'     => \$help
	     );
die $usage if $help;


sub parse_group_field {
    my( $group_field ) = @_;
    
    my ($igi, $gene_name, $native_id, $transcript_id, $exon_num, $exon_id);

    # Parse the group field
    foreach my $tag_val (split /;/, $group_field) {

        # Trim trailing and leading spaces
        $tag_val =~ s/^\s+|\s+$//g;

        my($tag, $value) = split /\s+/, $tag_val, 2;

        # Remove quotes from the value
        $value =~ s/^"|"$//g;
        $tag = lc $tag;

        if ($tag eq 'igi_id') {
            $igi = $value;
        }
        elsif ($tag eq 'gene_name') {
            $gene_name = $value;
        }
        elsif ($tag eq 'gene_id') {
            $native_id = $value;
        }
        elsif ($tag eq 'transcript_id') {
            $transcript_id = $value;
        }
        elsif ($tag eq 'exon_number') {
            $exon_num = $value;
        }
        elsif ($tag eq 'exon_id') {
            $exon_id = $value;
        }
        else {
            #warn "Ignoring group field element: '$tag_val'\n";
        }
    }
    return($igi, $gene_name, $native_id, $transcript_id, $exon_num, $exon_id);
}                                       # parse_group_field

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
    my ( @igis) = keys %$igi_hash;
    my $min = 1000000000000000000000;
    my $max = -1;
    my $sum = 0;
    my $minfeats = 1000000000000000000000;
    my $maxfeats = -1;
    my $sumfeats = 0;

    my $minexons = 1000000000000000000000;
    my $maxexons = -1;
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

my %igis_of_source;
my %native_ids_of_igi;

my @argv_copy = @ARGV; # may get gobbled up by the <> construct. 

### just read the complete file into the relevant hashes
GTF_LINE:
while (<>) {
    # taken largely from Bio::EnsEMBL::Utils::GTF_handler:
    next GTF_LINE if /^                 #/;
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
    my ($igi, $gene_name, $native_id, $transcript_id, $exon_num, $exon_id) =
      parse_group_field($group_field);
    
    unless ($igi) {
        warn("Skipping line with no igi_id: '$_'\n");
        next GTF_LINE;
    }
    unless ($native_id) {
        warn("Skipping line with no gene_id: '$_'\n");
        next GTF_LINE;
    }

    # keep track of marginals by looping over current one as well
    # as fictional source 'ALL':
    foreach my $s ($source, 'ALL') {
        #get previous record of this igi, if any:
        my ($nfeats, $min, $max, $nexons);
        if (defined $igis_of_source{$s}{$igi}) {
            ($nfeats, $min, $max, $nexons) = 
                @{$igis_of_source{$s}{$igi}};
        } else { 
            ($nfeats, $min, $max, $nexons)= (0, $start, $end, 0);
        }

        # keep track of the native gene_ids of a given igi in a given
        # source. As always, use a hash for faster collating (i.e. keys
        # %($native_ids_of_igi{$s}{$igi}) gives them all; the value is
        # irrelevant 
        $native_ids_of_igi{$s}{$igi}{$native_id}++;

        $min = $start if $start < $min;
        $max = $end if $end > $max;
        $nfeats++;
        $nexons = $exon_num if $exon_num > $nexons;

        # add record back in:
        $igis_of_source{$s}{$igi} = [$nfeats, $min, $max, $nexons];
        
        # pointless to keep track of exon statistics; call exon-lengths.awk
        # for that.
    }
}
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
    # print out native_id's clumped together in clusters of more than n
    foreach my $source ( 'ALL', @all_sources ) {
        print_chaining($source, $chaining);
    }
}

sub blurp {
    print '### $Id$  ', "\n";
    my (@stuff)  = ("on", `date`, "by", `whoami`, "@",  `hostname`);
    foreach (@stuff) {chomp};
    print "### run ", join(' ',@stuff), "\n";
    print "### argument(s): ", join(' ', @argv_copy), "\n";
    foreach (@argv_copy) {   print "### ", `ls -l $_`; }
    
    print "Sources: " , join( ' ', @all_sources), "\n";
}

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
            my @native_ids = keys %{$native_ids_of_igi{$source}{$igi}};
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
    blurp;

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
        #     my ($min, $max, $avg, $minfeats, $maxfeats, $avgfeats) = 
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
}                                       # all_stats

# print out native_id's clumped together in clusters of more than n
sub print_chaining {
    my ($source, $n)=@_;
    
    print "source $source:\n";
    for (my $i=$n; $i<int(@native_ids_per_igi_histo); $i++) {
        if ( $native_ids_per_igi_histo[ $i ]{$source}  ) { 
            print "    chaining together $i native ids:\n";
            my $igi_hash = $ {igis_per_chaining_group[ $i ]}{$source};
            foreach my $igi ( keys %$igi_hash) {
                print "      $igi: ";
                my @native_ids = keys %{$native_ids_of_igi{$source}{$igi}};
                print join(' ',@native_ids), "\n";
            }
        }
    }
}

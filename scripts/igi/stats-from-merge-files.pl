#!/usr/local/bin/perl

# $Id$ 
use strict;

sub parse_group_field {
    my( $group_field ) = @_;
    
    my ($igi, $gene_name, $gene_id, $transcript_id, $exon_num, $exon_id);

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
            $gene_id = $value;
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
    return($igi, $gene_name, $gene_id, $transcript_id, $exon_num, $exon_id);
}                                       # parse_group_field

sub gene_stats {
    my ($igi_hash) = @_;
    my ( @igis) = keys %$igi_hash;
    my $min = 1000000000000000000000;
    my $max = -1;
    my $sum = 0;
    my $minfeats = 1000000000000000000000;
    my $maxfeats = -1;
    my $sumfeats = 0;
    my $n =0; 
    foreach my $igi (@igis) {
        my ($nfeats, $start, $end) = @{$igi_hash->{$igi}};
        my $len = ($end - $start);
        $min = $len if $len <  $min ;
        $max = $len if $len >  $max ;
        $sum += $len;

        $minfeats = $nfeats if $nfeats <  $minfeats ;
        $maxfeats = $nfeats if $nfeats >  $maxfeats ;
        $sumfeats += $nfeats;

        $n++;

    }
    my ($avg, $avgfeats) = ('none', 'none');
    if ($n>0) {
        $avg = int($sum/$n);
        $avgfeats= int($sumfeats/$n);
    } 
    return ($min, $max, $avg, $minfeats, $maxfeats, $avgfeats);
}

my %igis_of_source;

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
    my ($igi, $gene_name, $gene_id, $transcript_id, $exon_num, $exon_id) =
      parse_group_field($group_field);
    
    unless ($igi) {
        warn("Skipping line with no igi_id: '$_'\n");
        next GTF_LINE;
    }

    # keep track of marginals by looping over current one and fictional
    # source 'all' at same time:
    foreach my $s ('ALL', $source) {
        #get previous record of this igi, if any:
        my ($nfeats, $min, $max);
        if (defined $igis_of_source{$s}{$igi}) {
            ($nfeats, $min, $max) = 
                @{$igis_of_source{$s}{$igi}};
        } else { 
            ($nfeats, $min, $max)= (0, $start, $end);
        }
# print "WAS: $start, $end, $min, $max\n";    
        $min = $start if $start < $min;
        $max = $end if $end > $max;
        $nfeats++;
# print "IS: $start, $end, $min, $max\n";    

        # add record back in:
        $igis_of_source{$s}{$igi} = [$nfeats, $min, $max];
    }
}

my @all_sources = sort keys %igis_of_source;
print "Found following sources:" , join( ' ', @all_sources), "\n";

print "number of igi's per source:\n";
foreach my $source (@all_sources) {
    my $igi_of_source = $igis_of_source{$source};
    my $num = int(keys %{$igi_of_source});
    print "$source: $num\n";
#    print "XXX ", join(':', keys %{$igi_of_source}), "\n";;
}
print "----\n";
# get rid of 'ALL':
@all_sources = (grep $_ ne 'ALL', sort keys %igis_of_source);

# compare igi's per source
print "pairwise overlaps:\n";
my %h = undef;
print "\t", join("\t\t", @all_sources) ,"\n";
SOURCE1:
foreach my $source1 (@all_sources) {
    my @igis_of_source1 = keys %{$igis_of_source{$source1}};
    my $source2;

    print "$source1\t";
  SOURCE2:
    foreach $source2 (@all_sources) {
        ## uncomment the undef for smaller output
        if ( $source2 gt $source1) { 
            print "-\t";
            next SOURCE2;
        }

        my ($shared, $tot) = (0,0); 
        my @igis_of_source2 = keys %{$igis_of_source{$source2}};
        foreach my $elt (@igis_of_source1) {
            if (defined $igis_of_source{$source2}{$elt}) {
                $shared++;
            }
            $tot++;
        }
        print "$shared/$tot\t";
    }
    print "\n"; 
}
print "----\n";
### sizes:
foreach my $source ('ALL', @all_sources) {
    my ($min, $max, $avg, $minfeats, $maxfeats, $avgfeats) = 
     gene_stats($igis_of_source{$source});

    print "gene stats on $source: ";
    print "minl=$min, maxl=$max, avgl=$avg, ";
    print "minf=$minfeats, maxf=$maxfeats, avgf=$avgfeats\n";
}
print "----\n";
## do overlaps
my @all_igis = keys %{$igis_of_source{'ALL'}};
my $n_igis = int(@all_igis);
my $n_sources = int(@all_sources);

# invert igis_of_source:
my %sources_of_igi= undef;
foreach my $igi (@all_igis) {
    foreach my $source (@all_sources ) { # excluding 'ALL', of course
        # does this source contain this igi?
        if (defined ${$igis_of_source{$source}}{$igi} ) {
          $sources_of_igi{$igi}{$source}++; # meaning: seen it
      }
    }
}

my @igis_of_n_sources = undef;

## do histogram: give numbers
foreach my $igi (@all_igis) {
    my $n  = int(keys %{$sources_of_igi{$igi}});
    my $loc = ${$igis_of_source{'ALL'}}{$igi}; # start,end etc. of this igi

#     ${$igis_of_n_sources[$n]}{$igi} = $loc this is for those that are in
# exactly two groups, but you want to know the  cumulation: simply add it 
# to all the clusterings. 
    for (my $i=$n; $i>=0; $i--) {
        ${$igis_of_n_sources[$i]}{$igi} = $loc
    }
}

# warn @igis_of_n_sources;
print "overlap totals (histogram)\n" ;
for(my $i = 0; $i<=$n_sources; $i++) {
    print "numbers of igis in at least $i sources: ",
      int(keys %{$igis_of_n_sources[$i]}), "/$n_igis\n";
}

print "----\n";
print "gene stats per cluster group:\n";
for(my $i = 1; $i<=$n_sources; $i++) {
    my ($min, $max, $avg, $minfeats, $maxfeats, $avgfeats) = 
        gene_stats($igis_of_n_sources[$i]) ;

    print "those in at least $i sources: ";
    print "minl=$min, maxl=$max, avgl=$avg, ";
    print "minf=$minfeats, maxf=$maxfeats, avgf=$avgfeats\n";
}



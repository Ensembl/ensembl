#!/usr/local/bin/perl

use strict;
# $Id$ 

# yet another little ad hoc script use Bio::EnsEMBL::Utils::GTF_handler;
# my $gtfh=Bio::EnsEMBL::Utils::GTF_handler->new(); # just for some parse
# functions!

# sub parse_group_field {
#     my ($_) = @_;
#     
#     my @ids = split ';';
# 
#     my $hash = undef;
#     foreach my $field ( @ids ) { 
#         my ($name, $value) = split ' ', $field;
#         $value =~ s/\"//g;
#         $hash->{$name}=$value;
#     } 
#     $hash;
# }

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

#         warn "XX: $tag $value\n";
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
    
    $igis_of_source{$source}{$igi}++;
# warn "source $source has :", keys %{$igis_of_source{$source}};

    # (more later)
    
#    $sources_of_igi{$igi}->{$source}++;
#    # (for easy statistics?)
# 
    
#     unless ($gene_id) {
#         warn("Skipping line with no gene_id: '$_'\n");
#         next GTF_LINE;
#     }
#     unless ($transcript_id) {
#         warn("Skipping line with no transcript_id: '$_'\n");
#         next GTF_LINE;
#     }
}

my @all_sources = keys %igis_of_source;
print "Found following sources:" , join( ' ', @all_sources), "\n";

print "total numbers:\n";
foreach my $source (@all_sources) {
    my $igi_of_source = $igis_of_source{$source};
    my $num = int(keys %{$igi_of_source});
    print "source $source: $num igi's\n";
}

# compare igi's per source
print "2-overlaps numbers:\n";
my %h = undef;
foreach my $source1 (@all_sources) {
    my $igis_of_source1 = $igis_of_source{$source1};
# warn "source $source1 has :", keys %{$igis_of_source{$source1}};
# warn "source1 $source1 has :", keys %$igis_of_source1;

    my @igis_of_source1 = keys %{$igis_of_source1};
    my $source2;
    foreach $source2 (@all_sources) {
        next if ($source2 le $source1);

        my $igis_of_source2 = $igis_of_source{$source2};
        my @igis_of_source2 = keys %{$igis_of_source2};
        
        foreach my $elt (@igis_of_source1) {
            $h{$elt} ++;
        }
        
        foreach my $elt (@igis_of_source2) {
            $h{$elt} += 2;
        }
    }

    foreach my $key (keys %h ) {
        warn "$key: $h{$key}\n";
    }


    #### This is broken., argghg!:
    print "Only in $source1 :" , int(map( ($_ == 1), (keys %h))) , "\n";
    print "Only in $source2 :" , int(map( ($_ == 2), (keys %h))) , "\n";
    print "In both $source1 and $source2 :" , int(map( ($_== 3),  (keys %h))) , "\n";
}

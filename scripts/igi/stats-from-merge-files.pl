#!/usr/local/bin/perl

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
    
    my ($igi_id, $gene_name, $gene_id, $transcript_id, $exon_num, $exon_id);

    # Parse the group field
    foreach my $tag_val (split /;/, $group_field) {

        # Trim trailing and leading spaces
        $tag_val =~ s/^\s+|\s+$//g;

        my($tag, $value) = split /\s+/, $tag_val, 2;

        # Remove quotes from the value
        $value =~ s/^"|"$//g;
        $tag = lc $tag;

        if ($tag eq 'igi_id') {
            $igi_id = $value;
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
    
    return($gene_name, $gene_id, $transcript_id, $exon_num, $exon_id);
}                                       # parse_group_field


GTF_LINE:
while (<>) {
    # taken largely from Bio::EnsEMBL::Utils::GTF_handler:
    next GTF_LINE if /^                 #/;
      next GTF_LINE if /^\s*$/;
    chomp;
    
    my @fields = split "\t", $_;
    my ($seq_name,    $source, $feature,
        $start,  $end,    $score,
        $strand, $phase,  $group_field)  = @fields;
    $feature = lc $feature;
    
    unless ($group_field) {
        warn("no group field: skipping : '$_'\n");
        next GTF_LINE ;
    }
    
    # Extract the extra information from the final field of the GTF line.
    my ($igi_id, gene_name, $gene_id, $transcript_id, $exon_num, $exon_id) =
      parse_group_field($group_field);
    
    unless ($igi_id) {
        warn("Skipping line with no igi_id: '$_'\n");
        next GTF_LINE;
    }
    
    ${$igi_ids_of_source{$source}}->{$igi_id}++;
    # (more later)
    
    ${$sources_of_igi_id{igi_id}}->{$source}++;
    # (for easy statistics)
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

my @all_sources = keys %igi_ids_of_source;
warn join ', ', @all_sources;

print "total numbers:\n";
foreach my $source (@all_sources) {
    my $igi_per_source = $igi_ids_of_source{$source};
    my $num = int(keys %{%igi_per_source});
    print "source $source: $num igi's\n";
}

print "2-overlaps numbers:\n";
foreach my $source (@all_sources) {
    my @igis_per_source = keys @{$igi_ids_of_source{$source}};
    foreach my $source2 (@all_sources) {
        next if $source2 le $source;
        my @igis_per_source2 = keys @{$igi_ids_of_source{$source2}};
        
        my %h = undef;

        foreach $elt (@igis_per_source) {
            $h{$elt} ++;
        }
        
        foreach $elt (@igis_per_source2) {
            $h{$elt} += 2;
        }
        
        print "Only in $source :" . int(grep( $_== 1, keys %h)) . "\n";
        print "Only in $source2 :" . int(grep( $_== 2, keys %h)) . "\n";
        print "In both $source and $source2 :" . int(grep( $_== 3,  keys %h)) . "\n";
    }
}

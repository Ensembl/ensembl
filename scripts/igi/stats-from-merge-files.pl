#!/usr/local/bin/perl

# $Id$ 

# yet another little ad hoc script 
use Bio::EnsEMBL::Utils::GTF_handler;
my $gtfh=Bio::EnsEMBL::Utils::GTF_handler->new(); # just for some parse functions!

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
    
    $gtf{$igi_id} .= "$source" unless 
    
#     unless ($gene_id) {
#         warn("Skipping line with no gene_id: '$_'\n");
#         next GTF_LINE;
#     }
#     unless ($transcript_id) {
#         warn("Skipping line with no transcript_id: '$_'\n");
#         next GTF_LINE;
#     }
}

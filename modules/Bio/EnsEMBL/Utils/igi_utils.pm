# $Id$

# some utils for the igi work
# copyright EnsEMBL (http://www.ensembl.org)
# Written by Philip lijnzaad@ebi.ac.uk

package Bio::EnsEMBL::Utils::igi_utils;

### put all the igi's found in the summary into one big hash and return
### it. The keys are the igi_ids; the values are a (ref to the) list of
### native ids
sub  read_igis_from_summary {
    my ($IN)  = @_;                     
    my %bighash = undef;

    SUMMARY_LINE:
    while (<$IN>) {
        next SUMMARY_LINE if /^#/;
          next SUMMARY_LINE if /^\s*$/;
        chomp;
    
        my @fields = split "\t", $_;
        my ($seq_name, $source, $feature,
            $start,  $end,    $score,
            $strand, $phase,  $group_field)  = @fields;
        $feature = lc $feature;
        
        unless ($group_field) {
            warn("no group field: skipping : '$_'\n");
            next SUMMARY_LINE ;
        }
        
        # Extract the extra information from the final field of the GTF line.
        my ($igi, $gene_name, $native_ids, $transcript_id, $exon_num, $exon_id) =
          parse_group_field($group_field);

        unless ($igi) {
            warn("no  igi; skipping: '$_'\n");
        }
        
        $big_hash{$igi} = $native_ids;
    }                                   # while <$IN>
    return \%big_hash;
}                                       # read_igis_from_summary

### put all the transcript_ids found in the summary into one big hash and
### return it. Key is transcript id; value is the summed length of all its
### features.
sub read_transcript_ids_from_summary {
    my ($IN)  = @_;                     
    my %bighash = undef;

    SUMMARY_LINE:
    while (<$IN>) {
        next SUMMARY_LINE if /^#/;
          next SUMMARY_LINE if /^\s*$/;
        chomp;
    
        my @fields = split "\t", $_;
        my ($seq_name, $source, $feature,
            $start,  $end,    $score,
            $strand, $phase,  $group_field)  = @fields;
        $feature = lc $feature;
        
        unless ($group_field) {
            warn("no group field: skipping : '$_'\n");
            next SUMMARY_LINE ;
        }
        
        # Extract the extra information from the final field of the GTF line.
        my ($igi, $gene_name, $native_ids, $transcript_ids, $exon_num, $exon_id) =
          parse_group_field($group_field);

        if (@$transcript_ids==0) {
            die("no transcript_ids; skipping: '$_'\n");
        }

        if (@$transcript_ids==0) {
            die("no transcript_ids: '$_'\n");
        }

        foreach my $tid ( @$transcript_ids ) {
            $big_hash{$tid} ++;
        }
        
    }                                   # while <$IN>
    return \%big_hash;
}                                       # read_transcript_ids_from_summary


# return a list of fields from a gtf file. 
# native_id and trancript_id are lists, since the 
# summary files usually contain more than one native id and transcript_id
sub parse_group_field {
    my( $group_field ) = @_;
    
    my ($igi, $gene_name, @native_ids, @transcript_ids, $exon_num, $exon_id);

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
            push @native_ids,  $value;
        }
        elsif ($tag eq 'transcript_id') {
            push @transcript_ids, $value;
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
    return($igi, $gene_name, \@native_ids, \@transcript_ids, $exon_num, $exon_id);
}                                       # parse_group_field


### sorts igi3_5 igi3_10 on the basis of the numbers.  call as :
###
###    sort { Bio::EnsEMBL::Utils::igi_utils::by_igi_number($a, $b) } @igis; 
### 
### 
sub by_igi_number {
  my ($aa, $bb) = (substr($_[0],5), substr($_[1],5));
  $aa <=> $bb;
}

1;

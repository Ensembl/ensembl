#
#
# Cared for by Ensembl  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::KnownMapping::Pmatch

=head1 SYNOPSIS


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::KnownMapping::Pmatch;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::FeatureFactory;
use Bio::EnsEMBL::KnownMapping::MappingConf qw (
						QUERY
						KNOWN_PROTEINS_FASTA
						PMATCH_OUT
						TARGET_IDT
						QUERY_IDT
						);



sub new {
    my $class= shift;
    my $self = {};
    bless $self,$class;
   
    # make sure that's the basic options are defined
    $self->throw("One or few of the following options has not been defined:\nQUERY: $QUERY\nKNOWN_PROTEINS_FASTA: $KNOWN_PROTEINS_FASTA\nTARGET_IDT: $TARGET_IDT\nQUERY_IDT: $QUERY_IDT\n") 
    unless (($QUERY) && ($KNOWN_PROTEINS_FASTA) && ($TARGET_IDT) && ($QUERY_IDT));

    return $self;
}


sub run {
    my $self = shift;
    print STDERR "Running pmatch\n";

    my $runpmatch = "pmatch -T 14 $KNOWN_PROTEINS_FASTA $QUERY > $$.pmatch";
    system "$runpmatch";

}

sub parse {
    my $self = shift;
    my %raw_features;

#pmatch raw output is as follow
#length, qid, qstart, qend, qperc, tid, tstart, tend, tperc

    print STDERR "Sorting pmatch\n";
    
    open (IN,"$$.pmatch") || die "Can't open pmatch intermediary file: $$.pmatch\n";
    
    while (<IN>) {
	#print "$_";
	chomp;
	
	my (@array) = split;

	my $length = $array[0];
	my $id = $array[1];
	my $hid = $array[5];
	my $start = $array[2];
	my $end = $array[3];
	my $hstart = $array[6];
	my $hend = $array[7];
	my $perc = $array[4];
	my $hperc = $array[8];

	my $new_id = $hid.":".$id;
	

#	 my ($self,$start,$end,$strand,$score,$source,$primary,$seqname,$hstart,$hend,
#        $hstrand,$hscore, $hsource,$hprimary,$hseqname, $e_value, $perc_id, 
#        $phase, $end_phase) = @_;

	my $fp      = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
	
	$fp->set_all_fields($start,$end,1,$perc,"NULL",$new_id,$id,$hstart,$hend,1,$hperc,"NULL","NULL",$hid,$length);

	push (@{$raw_features{$new_id}},$fp);
    }

    foreach my $k (keys %raw_features) {
	my $score;
	my $idt;
	my $hidt;
	
	my @array = @{$raw_features{$k}};
	

	foreach my $a (@array) {
	    print $a->p_value."\t".$a->start."\t".$a->end."\t".$a->id."\t".$a->hstart."\t".$a->hend."\t".$a->feature2->id."\n";
	}
	
	print "BEF: ".scalar(@array)."\n";

	my @features_filtered = $self->_greedy_filter(@array);

	@features_filtered = sort {$a->hstart <=> $b->hstart} @features_filtered;



#my ($len,$qid,$qstart,$qend,$qperc,$qlen,$tid,$tstart,$tend,$tperc,$tlen) = split /\t/;	

	open (INT,">$$.intermediate") || die "Can't open intermediate\n";

	foreach my $b (@features_filtered) {
	    print INT $b->p_value."\t".$b->start."\t".$b->end."\t".$b->id."\t".$b->hstart."\t".$b->hend."\t".$b->feature2->id."\n";
	}

	print "AFT: ".scalar(@features_filtered)."\n\n";
	
    }

    $self->process_matches(*IN); 

}


sub _greedy_filter {
    my ($self,@features) = @_;

    @features = sort {$b->p_value <=> $a->p_value} @features;
    
    my @features_tmp;
    
    foreach my $fp (@features) {
	if (! scalar @features_tmp) {
	    push @features_tmp, $fp;
	    next;
	}
	
	my $add_fp = 1;
	foreach my $feature_tmp (@features_tmp) {
	    my ($start,$end,$hstart,$hend) = ($feature_tmp->start,$feature_tmp->end,$feature_tmp->hstart,$feature_tmp->hend);


##############################################################################	    
#This step allows to get rid of repeated sequence which match everywhere, eg:#
#      24      38              654     668                                   #
#      24      39              653     668                                   #
#      25      40              653     668                                   #
#      26      41              653     668                                   #
#      24      37              655     668                                   # 
#      27      42              653     668                                   #
#      28      43              653     668                                   #
#      29      44              653     668                                   #
##############################################################################


	    unless ((($fp->hstart < $hstart) || ($fp->hend > $hend)) && (($fp->start < $start) || ($fp->end > $end)))   {
		$add_fp = 0;
	    }

	   

	}
	push @features_tmp, $fp if ($add_fp);
    }

    
    #@features_tmp = sort {$a->hstart <=> $b->hend} @features_tmp;


    my @features_filtered;

    
    foreach my $ff (@features_tmp) {
	if (! scalar @features_filtered) {
	    push @features_filtered, $ff;
	    next;
	}


	foreach my $ft (@features_filtered) {
	    my ($start,$end,$hstart,$hend) = ($ft->start,$ft->end,$ft->hstart,$ft->hend);
	    
#This check if the end of the feature does not overlap with other feature. First check for hend. If there is any overlap take out the difference, do the same for end.
	    if ($ff->hend >= $hstart && $ff->hend <= $hend) {
		my $diff = $ff->hend - $hstart + 1;
		$ff->hend($hstart-1);
		$ff->end($ff->end - $diff);
		
		if ($ff->end >= $start && $ff->end <= $end) {
		    my $diff1 = $ff->end - $start + 1;
		    $ff->end($start-1);
		    $ff->hend($ff->hend - $diff);
		}
		
	    }

#Do the same but check for hstart	    
	      if ($ff->hstart >= $hstart && $ff->hstart <= $hend) {
		my $diff = $hend - $ff->start + 1;
		$ff->hstart($hstart+1);
		$ff->start($ff->start + $diff);
		
		if ($ff->start >= $start && $ff->start <= $end) {
		    my $diff1 = $hend - $ff->start + 1;
		    $ff->start($start+1);
		    $ff->hstart($ff->hstart + $diff);
		}
		
	    }
	    
	}

	push @features_filtered,$ff;
	
    }
	return @features_filtered;
} 


sub _process_matches {
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
        # get the total matching sequence (pmatch gives multiple hits either due to e.g.
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
	
	print  OUT $top->name."\t$know\tcandidate\t".$top->qperc."\t".$top->tperc."\n";

	foreach my $repeat (@array) {
	    if( $repeat->tperc+1 >= $top->tperc ) {
		print  OUT $repeat->name."\t$know\tduplicate\t".$repeat->qperc."\t".$repeat->tperc."\n";
	    }
	}
    }
}



1;





#
# BioPerl module for Bio::EnsEMBL::Analysis::CdnaResolver
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::CdnaResolver - Resolves cdnas onto genomic DNA

=head1 SYNOPSIS
 
    $resolver = Bio::EnsEMBL::Analysis::CdnaResolver->new();
    $resolver->log($log); # log is Bio::EnsEMBL::Analysis::Log object

    @seqfeature = $resolver->resolve($annseq,@cdna_accessions); # Accession numbers of cdnas

    # seqfeatures are a list of new seqfeatures which can be optionally added to
    # to annseq

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Analysis::CdnaResolver;
use vars qw($AUTOLOAD @ISA);
use strict;

use Bio::AnnSeq;
use Bio::AnnSeqIO;

use Bio::Tools::Sim4::Results;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 resolve

 Title   : resolve
 Usage   : $resolver->resolve($aseq,@cdna_accessions)
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub resolve{
    my ($self,$aseq,@cdna_accessions) =	@_;
    my @exon_set_out;

    # write out the sequence somewhere sensible.
    my %top_sf;

    my $nseq = $aseq->seq->seq();
    my $bnseq = Bio::Seq->new( -seq => lc($nseq), -type => 'Dna', -id => $aseq->seq->id() );
    
    my $seqout1 = Bio::SeqIO->new( -format => 'Fasta', -file => ">/tmp/aseq.$$");
    $seqout1->write_seq($bnseq);
    $seqout1 = undef;
    
    $bnseq = $nseq = undef;

  ACC: 
    foreach my $embl_acc ( @cdna_accessions ) {
	
	my $cdna = $self->_fetch_embl_acc($embl_acc);
	
	# skip if we can't load it or no cdna.
	if( !defined $cdna) {
	    $self->_message("Unable to read $embl_acc. Skipping. (Apologies)");
	    next ACC;
	}

	# get out CDS coordinates.
	
	my $cds_start;
	my $cds_end;
	
	foreach my $sf ( $cdna->top_SeqFeatures ) {
	    if( $sf->primary_tag() eq 'CDS' ) {
		$cds_start = $sf->start;
		$cds_end   = $sf->end;
	    }
	}
	
	if( !defined $cds_start ) {
	    $self->_message("EMBL $embl_acc has not single CDS line. Skipping");
	}
	
    
	# wrap everything in an eval to trap exceptions
	eval {
	    # assumme file is in fasta format for sim4 stuff.
	    
	    my $seqout = Bio::SeqIO->new( '-format' => 'Fasta', -file => ">/tmp/sim4.cdna.$$");
	    $cdna->seq->setseq(lc($cdna->seq->seq));
	    $seqout->write_seq($cdna->seq);
	    $seqout = undef;
	    print STDERR "file /tmp/aseq.$$ vs /tmp/sim4.cdna.$$ \n";

	    my $sim4 = Bio::Tools::Sim4::Results->new( -file => "sim4 /tmp/aseq.$$ /tmp/sim4.cdna.$$ |" );
	    

	    foreach my $es ( $sim4->each_ExonSet ) {
		
		# save this for returning to ExonSets.
		push(@exon_set_out,$es);
		
		# we need to process this into a set of CDS regions. 
		# we need to get the first exon.
		
		my $seen_start = 0;
		my @exons = $es->each_Exon;

		if( $#exons == -1 ) {
		    die("No exons in $embl_acc sim4 run!");
		}

		
		# sort by whether homology is start/end.
		
		@exons = sort { $a->homol_SeqFeature->start <=> $b->homol_SeqFeature->end } @exons;
		
		my $strand;
		
		$strand = $exons[0]->homol_SeqFeature->strand();
		
		my $exon; # temporay variable of the current exon
		while( $exon = shift @exons ) {
		    if( $exon->homol_SeqFeature->start < $cds_start && 
			$exon->homol_SeqFeature->end > $cds_start ) {
			$seen_start =1;
			last;
		    }
		}

		if( $seen_start == 0 ) {
		    $self->_message("Have not seen start. Skipping an exon set");
		    next;
		}
		
	    
		#
		# Build a top level Generic SeqFeature object
		# and add in the individual exons
		#
	    
		my $top_sf = Bio::SeqFeature::Generic->new();
		$top_sf->primary_tag('CDS_span');
		$top_sf->source_tag('Sim4_CDS');
		$top_sf->_parse->{'parent_homogenous'} =1;
		$top_sf->strand($strand);
		$top_sf->add_tag_value('match_cdna',$embl_acc);

		# first exon. Has to change start.
		my $sf = Bio::SeqFeature::Generic->new();
		if( $strand == 1 ) {
		    $sf->start($exon->start + ($cds_start - $exon->homol_SeqFeature->start));
		    $sf->end($exon->end);
		} else {
		    $sf->start($exon->start);
		    $sf->end($exon->end - ($cds_start - $exon->homol_SeqFeature->start));
		}
		
		$sf->strand($strand);
		$sf->primary_tag('CDS');
		$sf->source_tag('Sim4_CDS');
		$sf->add_tag_value('match_cdna',$embl_acc);
		$top_sf->add_sub_SeqFeature($sf,'EXPAND');
		
		my $seen_end = 0;
		my $end_exon;
		foreach $exon ( @exons ) {
		    # the remainder
		    if( $exon->homol_SeqFeature->start() < $cds_end && 
			$exon->homol_SeqFeature->end() > $cds_end ) {
			$seen_end = 1;
			$end_exon = $exon;
			last;
		    }

		    my $sf = Bio::SeqFeature::Generic->new();
		    if( $strand == 1 ) {
			$sf->start($exon->start);
			$sf->end($exon->end);
		    } else {
			$sf->start($exon->start);
			$sf->end($exon->end);
		    }

		    $sf->strand($strand);
		    $sf->primary_tag('CDS');
		    $sf->source_tag('Sim4_CDS');
		    $sf->add_tag_value('match_cdna',$embl_acc);
		    $top_sf->add_sub_SeqFeature($sf,'EXPAND');
		}

		$exon= $end_exon;

		# last exon.
		if( $seen_end == 0 ) {
		    $self->_message("Did not see end exon. Yuk.\n");
		} else {
		    
		    my $sf = Bio::SeqFeature::Generic->new();
		    if( $strand == 1 ) {
			$sf->start($exon->start);
			$sf->end($exon->start + ($cds_end - $exon->homol_SeqFeature->start));
		    } else {
			$sf->start($exon->end - ($cds_end - $exon->homol_SeqFeature->start));
			$sf->end($exon->end);
		    }
		    
		    $sf->strand($strand);
		    $sf->primary_tag('CDS');
		    $sf->source_tag('Sim4_CDS');
		    $sf->add_tag_value('match_cdna',$embl_acc);
		    $top_sf->add_sub_SeqFeature($sf,'EXPAND');
		}

		# add in this feature
		
		$top_sf->strand($strand);
		$top_sf{$top_sf} = $top_sf;
		
	    }
	};

	# catch exceptions
	
	if( $@ ) {
	    $self->_message("unable to process $embl_acc\n$@\n");
	}
	
    }

    #
    # Sort by size
    #

    my @size = sort { $a->length <=> $b->length } values %top_sf;

    #
    # Go down the list, marking off guys which are subsummed by deleting from the hash
    #
    
    foreach my $sf ( @size ) {
	foreach my $test ( values %top_sf ) {
	    if( $test == $sf ) {
		next;
	    }
	    
	    if( compare_top_seqfeature($test,$sf) == 1 ) {
		# remove sf
		delete $top_sf{$sf};
	    }
	}
    }

    return (@exon_set_out,values %top_sf);

    
}

=head2 fetch_annseq_func

 Title   : fetch_annseq_func
 Usage   : $obj->fetch_annseq_func($newval)
 Function: 
 Returns : value of fetch_annseq_func
 Args    : newvalue (optional)


=cut

sub fetch_annseq_func{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'fetch_annseq_func'} = $value;
    }
    return $obj->{'fetch_annseq_func'};

}

=head2 _fetch_embl_acc

 Title   : _fetch_embl_acc
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _fetch_embl_acc{
   my ($self,$acc) = @_;

   if( $self->fetch_annseq_func ) {
       return &{$self->fetch_annseq_func}($acc);
   } else {
       	my $aseqio = Bio::AnnSeqIO->new( '-format' => 'EMBL', 
					 -file => "efetch -a em:$acc |");    
	my $cdna;
	eval {
	    $cdna = $aseqio->next_annseq();
	};
	if( $@ ) {
	    $self->_message("unable to fetch annseq object for $acc!\n$@\n");
	    return undef;
	} else {
	    return $cdna;
	}
    }
}

=head2 _message

 Title   : _message
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _message{
   my ($self,$mess) = @_;
   
   print STDERR "Resolver: $mess\n";


}

##################
# SubRoutines
##################

sub compare_top_seqfeature {
    my ($one,$two) = @_;

    # is two completely contained by one

    if( $one->strand ne $two->strand ) {
	return -1;
    }
    my $strand = $one->strand;

    if( $one->start > $two->start || $one->end < $two->end ) {
	# can't be contained
	return -1;
    }
    
    # first and last sub features handled differently.

    my @exons = $two->sub_SeqFeature;
    my $first = shift @exons;
    my $last  = pop @exons;

    my @matched_exons = $one->sub_SeqFeature;

    #
    # Order is important. It is not good enough just that
    # exons in two match one, as if one has an additional exon,
    # then they are separate
    #

    while ( my $current_exon = shift @matched_exons ) {
	if( $strand == 1 ) {
	    if( $current_exon->start < $first->start && $current_exon->end > $first->start ) {
		if( $first->end != $current_exon->end ) {
		    return -1;
		} else {
		    last;
		}
	    }
	} else {
	    if( $current_exon->end > $first->end && $current_exon->start < $first->end ) {
		if( $first->start != $current_exon->start ) {
		    return -1;
		} else {
		    last;
		}
	    }
	}
    }

    if( $#matched_exons == -1 ) {
	return -1;
    }

    #
    # ok. Now each exon should agree 
    #

    while( my $current_two_exon = shift @exons ) {
	my $current_exon = shift @matched_exons;
	if( !defined $current_exon ) {
	    return -1;
	}

	if( $current_two_exon->start != $current_exon->start ||
	    $current_two_exon->end != $current_exon->end ) {
	    return -1;
	}
    }

    #
    # final exon
    #

    my $last_one_exon = shift @matched_exons;

    if( $strand == 1 ) {
	if( $last->start != $last_one_exon->start || $last->end > $last_one_exon->end )  {
	    return -1;
	}
    } else {
	if( $last->end != $last_one_exon->end || $last->start < $last_one_exon->start )  {
	    return -1;
	}
    }

    return 1;
}

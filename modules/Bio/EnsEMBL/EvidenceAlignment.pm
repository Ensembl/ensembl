#
#
# Cared for by Ensembl <ensembl-dev@ebi.ac.uk>
#
# You may distribute the source code of the module
# under the same terms as perl itself.
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::EvidenceAlignment.pm

=head1 SYNOPSIS

 my $ea = Bio::EnsEMBL::EvidenceAlignment->new(
                          -DBADAPTOR    => $dba,
                          -TRANSCRIPTID => $tr_stable_id);
 my $seqs_arr_ref = $ea->fetch_alignment;
 $ea->transcriptid($other_tr_stable_id);
 my $other_seqs_arr_ref = $ea->fetch_alignment;

=head1 DESCRIPTION

Gives a transcript and its evidence as an alignment, with padding with
"-" as required. Proteins are given before DNA. Coordinates in the
database are used to recreate the alignment, and must be correct or some
data may not be displayed.

=head1 CONTACT

Ensembl: ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::EvidenceAlignment;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::SeqIO;
use Bio::EnsEMBL::Gene;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

=head2 new

    Title   :   new
    Usage   :   my $ea = Bio::EnsEMBL::EvidenceAlignment->new(
                                         -DBADAPTOR    => $dba,
                                         -TRANSCRIPTID => $tr_stable_id);

    Function:   Initialises EvidenceAlignment object
    Returns :   An EvidenceAlignment object
    Args    :   Database adaptor object and transcript stable ID string

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($transcriptid, $dbadaptor) = $self->_rearrange(['TRANSCRIPTID',
                                                      'DBADAPTOR'],
						      @args);
  $self->transcriptid($transcriptid);
  $self->dbadaptor($dbadaptor);

  return $self; # success - we hope!
}

=head2 dbadaptor

    Title   :   dbadaptor
    Usage   :   $ea->dbadaptor($dba);
    Function:   get/set for database adaptor object

=cut

sub dbadaptor {
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{'db_adaptor'} = $value;
  }
  return $obj->{'db_adaptor'};
}

=head2 transcriptid

    Title   :   transcriptid
    Usage   :   $ea->transcriptid($transcript_stable_id);
    Function:   get/set for transcript stable id string

=cut

sub transcriptid {
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{'transcript_id'} = $value;
  }
  return $obj->{'transcript_id'};
}

=head2 fetch_alignment

    Title   :   fetch_alignment
    Usage   :   my $seq_arr_ref = $ea->fetch_alignment;
    Function:   gets aligned transcript and evidence
    Returns :   reference to array of Bio::PrimarySeq

=cut

sub fetch_alignment {
  my ($self) = @_;
  $self->throw("must have a transcript stable ID and a DB adaptor object")
    unless ($self->transcriptid && $self->dbadaptor);
  return $self->_get_aligned_evidence($self->transcriptid, $self->dbadaptor);
}

# get a sequence from cache if we have it, otherwise from seqfetcher
sub _get_seq {
  my ($self, $seqfetcher, $cache_arr_ref, $name) = @_;
  $self->throw("interface fault") if (@_ != 4);

  foreach my $cache_entry (@$cache_arr_ref) {
    if ($$cache_entry{'name'} eq $name) {
      return $$cache_entry{'seqobj'};
    }
  }
  my %new_entry = ();
  $new_entry{'seqobj'} = $seqfetcher->run_pfetch($name);
  $new_entry{'name'} = $name;
  push @$cache_arr_ref, \%new_entry;
  return $new_entry{'seqobj'};
}

sub _get_transcript_nuc {
  my ($self, $exon_arr_ref) = @_;
  $self->throw("interface fault") if (@_ != 2);

  my $retval = "";
  my $seq_str;
  for (my $i = 0; $i <= $#$exon_arr_ref; $i++) {
    my $exon_seq = $$exon_arr_ref[$i]->seq->seq;
    if ($i & 1) {
      $seq_str = "\L$exon_seq";
    } else {
      $seq_str = "\U$exon_seq";
    }
    $retval .= $seq_str;
  }
  return $retval;
}

sub _get_transcript_pep {
  my ($self, $exon_pep_arr_ref) = @_;
  $self->throw("interface fault") if (@_ != 2);
  
  my $retval = "";
  my $seq_str;
  for (my $i = 0; $i <= $#$exon_pep_arr_ref; $i++) {
    my $pep = $$exon_pep_arr_ref[$i];
    if ($i & 1) {
      $seq_str = "\L$pep";
    } else {
      $seq_str = "\U$pep";
    }
    $retval .= $seq_str;
  }
  return $retval;
}

# _get_aligned_evidence: takes a transcript ID and a DB adaptor
# returns ref to an array of Bio::PrimarySeq

sub _get_aligned_evidence {
  my ($self, $transcript_id, $db) = @_;
  $self->throw("interface fault") if (@_ != 3);

  my $sgp = $db->get_StaticGoldenPathAdaptor;
  my @evidence_arr;	# a reference to this is returned
  my $evidence_obj;

  my $ta = $db->get_TranscriptAdaptor;
  my $ga = $db->get_GeneAdaptor;
  my $ea = $db->get_ExonAdaptor;
  my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher->new();
  $seqfetcher->pfetch("/usr/local/pubseq/bin/pfetch");

  # get all exons off a VC
  my $gene = $ga->fetch_by_transcript_stable_id($transcript_id);
  my $vc = $sgp->fetch_VirtualContig_of_gene($gene->stable_id, 100);
  my @genes_in_vc = $vc->get_Genes_by_Type('ensembl');
  my $transcript_obj;
  GENE_LOOP:
  foreach my $gene_in_vc (@genes_in_vc) {
    if ($gene_in_vc->stable_id eq $gene->stable_id)
    {
      my @transcripts_in_vc = $gene_in_vc->each_Transcript;
      foreach my $transcript_in_vc (@transcripts_in_vc) {
        if ($transcript_in_vc->stable_id eq $transcript_id) {
          $transcript_obj = $transcript_in_vc;
          last GENE_LOOP;
        }
      }
    }
  }

  # get transcript start and end on VC
  my @all_exons = $transcript_obj->get_all_Exons;
  my $transcript_start = 100000000;	# v. large value
  my $transcript_end = 0;
  foreach my $exon (@all_exons) {
    if ($exon->start < $transcript_start) {
      $transcript_start = $exon->start;
    }
    if ($exon->end > $transcript_end) {
      $transcript_end = $exon->end;
    }
  }

  # protein evidence

  my @utr_free_exons = $transcript_obj->translateable_exons;
  my $translation = $transcript_obj->translate->seq;
  print STDERR "XXX translation length ".length($translation),"\n";
  my $nucseq = $self->_get_transcript_nuc(\@all_exons);
  my $nuc_tran_len = length $nucseq;

  # record whether each exon is totally UTR, simply a frameshift (hence
  # not to be translated), or at least partly translating
  my @exon_status = ();
  for (my $i = 0; $i <= $#all_exons; $i++) {
    print STDERR "XXX exon $i len ".($all_exons[$i]->end - $all_exons[$i]->start + 1) . "\n";
    $exon_status[$i] = 'utr';
  }
  for (my $i = 0; $i <= $#all_exons; $i++) {
    if (($all_exons[$i]->end - $all_exons[$i]->start + 1) < 3) {
      $exon_status[$i] = 'frameshift';
    } else {
      foreach my $utr_free_exon (@utr_free_exons) {
        if (($utr_free_exon->start >= $all_exons[$i]->start) 
        && ($utr_free_exon->end <= $all_exons[$i]->end)) {
          $exon_status[$i] = 'translating';
        }
      }
    }
  }

  # get length of the UTR witin the first exon (excludes length of any
  # exons that fall entirely within the 5' UTR)
  # and also get the total 5' UTR length
  my $total_utr_len = 0;
  my $utr_region_in_first_exon = 0;
  UTR_OUTER:
  for (my $i = 0; $i <= $#all_exons; $i++) {
    if ($exon_status[$i] eq 'translating') {
      foreach my $utr_free_exon (@utr_free_exons) {
        if ($utr_free_exon->overlaps($all_exons[$i])) {
          if ($all_exons[$i]->strand > 0) {	# fwd gene
  	    if ($utr_free_exon->start != $all_exons[$i]->start) {
	      $utr_region_in_first_exon = $utr_free_exon->start
	                                  - $all_exons[$i]->start;
	      last UTR_OUTER;
	    }
	  } else {	# rev gene
	    if ($utr_free_exon->end != $all_exons[$i]->end) {
	      $utr_region_in_first_exon = $all_exons[$i]->end
	                                  - $utr_free_exon->end;
	      last UTR_OUTER;
	    }
	  }
        }
      }
    } else {
      $total_utr_len += $all_exons[$i]->end - $all_exons[$i]->start + 1;
    }
  }
  $total_utr_len += $utr_region_in_first_exon;
  print STDERR "XXX total_utr_len $total_utr_len\n";
  print STDERR "XXX utr_region_in_first_exon $utr_region_in_first_exon\n";

  # create an array of all exons that at least partly translate
  print STDERR "XXX exon cnt ".($#all_exons+1)."\n";
  my @exons_to_display = ();
  for (my $i = 0; $i <= $#all_exons; $i++) {
    if ($exon_status[$i] eq 'translating') {
      push @exons_to_display, $all_exons[$i];
    } else {
      print STDERR "XXX exon $i not translating\n";
    }
  }

  my $evidence_line = $translation;
  my $prot_tran_len = length $translation;
  # pad between residues, for comparability with codons
  my $space_free_length = length $evidence_line;
    my $padded_evidence_line = "";
    for (my $i = 0; $i < $space_free_length; $i++) {
      $padded_evidence_line .= substr $evidence_line, $i, 1;
      $padded_evidence_line .= "--";
    }
    # adjust for 3' and 5' UTRs
    $evidence_line = ('-' x $total_utr_len) . $padded_evidence_line;
    while (length($evidence_line) < $nuc_tran_len) {
      $evidence_line .= '-';
    }

  $evidence_obj = Bio::PrimarySeq->new(
                    -seq              => $evidence_line,
                    -id               => 0,
     		    -accession_number => $transcript_obj->stable_id,
		    -moltype          => 'protein'
		  );

  push @evidence_arr, $evidence_obj;

  my @seqcache = ();

  #my $total_exon_len = $utr_region_in_first_exon;
  my $total_exon_len = 0;
  my @pep_evidence_arr = ();
  for (my $i = 0; $i <= $#exons_to_display; $i++) {
    my $start = $exons_to_display[$i]->start;
    my @features = $exons_to_display[$i]->each_Supporting_Feature;

    my $last_feat = undef;
    PEP_FEATURE_LOOP:
    foreach my $feature(@features) {
      next PEP_FEATURE_LOOP if ($last_feat
      && ($last_feat->start == $feature->start)
      && ($last_feat->end == $feature->end)
      && ($last_feat->hstart == $feature->hstart)
      && ($last_feat->hseqname eq $feature->hseqname));
      my $fstart = $feature->start;
      my $hit_seq_obj = $self->_get_seq($seqfetcher, \@seqcache,
                                        $feature->hseqname);
      next PEP_FEATURE_LOOP if (! $hit_seq_obj);
      if ($hit_seq_obj->moltype eq "protein") {
        my $hlength = $feature->hend - $feature->hstart + 1;
	print "protein hit length (aa residues): $hlength\n";
        my $length = $feature->end - $feature->start + 1;
	print "protein feature len (bases)     : $length\n";
        next PEP_FEATURE_LOOP unless ($length == (3 * $hlength));
	my $hseq;
	eval {
          $hseq = substr $hit_seq_obj->seq, $feature->hstart - 1, $hlength;
	};
	next PEP_FEATURE_LOOP if ($@);	# data wrong
        my $hindent;
        if ($exons_to_display[$i]->strand > 0) {
          $hindent = ($total_exon_len + $fstart - $start) / 3;
        } else {
          $hindent = ($total_exon_len + $exons_to_display[$i]->end
	             - $feature->end) / 3;
        }
        $hindent = int $hindent;
        if ($i == 0) {
          $hindent -= $utr_region_in_first_exon;
        }
	if ($hindent < 0) {
	  $hindent = 0;
	}
        my %hit_details = ( 'hseqname'    => $feature->hseqname,
                            'hstart'      => $feature->hstart,
		            'hend'        => $feature->hend,
			    'hlength'     => $hlength,
			    'hseq'        => $hseq,
			    'hindent'     => $hindent,
			    'score'       => $feature->score,
			    'start'       => $start,
			    'end'         => $feature->end,
			    'exon_strand' => $exons_to_display[$i]->strand,
			    'exon_length' => $exons_to_display[$i]->end
			                     - $start + 1,
  			    'exon'        => $exons_to_display[$i]
			  );
        push @pep_evidence_arr, \%hit_details;
      }
      $last_feat = $feature;
    }
    $total_exon_len += $exons_to_display[$i]->end
    - $exons_to_display[$i]->start + 1; # 3 * length($exon_peps[$i]);
  }

  my @sorted_pep_evidence_arr = sort {   ( $$a{'hseqname'} =~ /^ENST/
                                           ? -1 : 0 )
			              || ( $$b{'hseqname'} =~ /^ENST/
                                           ? 1 : 0 )
				      || ( $$a{'hseqname'} cmp $$b{'hseqname'} )
				      || ( $$a{'hindent'} <=> $$b{'hindent'} )
 			        } @pep_evidence_arr;

  $evidence_line = "";
  my $uppercase = 1;	# case of sequence for output
  my $hit = $sorted_pep_evidence_arr[0];
  my $prev_exon = $$hit{'exon'};
  my $prev_hseqname = "-" x 1000;	# fake initial ID
  for (my $i = 0; $i <= $#sorted_pep_evidence_arr; $i++) {
    my $hit = $sorted_pep_evidence_arr[$i];
    if ($$hit{'exon'} ne $prev_exon) {
      $uppercase = ! $uppercase;
    }
    my $seq_str = $$hit{'hseq'};
    if ($uppercase) {
      $seq_str = "\U$seq_str";
    } else {
      $seq_str = "\U$seq_str";
      # $seq_str = "\L$seq_str";
    }
    
    if ($$hit{'hseqname'} ne $prev_hseqname) {	# make new evidence line
      $evidence_line = "-" x $prot_tran_len;
    }

    # splice in the evidence fragment
    substr $evidence_line, $$hit{hindent}, length $seq_str, $seq_str;

    # store if end of evidence line

    if (($i == $#sorted_pep_evidence_arr) ||
        ($sorted_pep_evidence_arr[$i+1]{'hseqname'} ne $$hit{'hseqname'}))
    {
      $evidence_line = substr $evidence_line, 0, $prot_tran_len;

      # pad between residues, for comparability with codons
      my $space_free_length = length $evidence_line;
      my $padded_evidence_line = "";
      for (my $i = 0; $i < $space_free_length; $i++) {
        $padded_evidence_line .= substr $evidence_line, $i, 1;
	$padded_evidence_line .= "--";
      }

      # adjust for 3' and 5' UTRs
      $evidence_line = ('-' x $total_utr_len) . $padded_evidence_line;
      while (length($evidence_line) < $nuc_tran_len) {
        $evidence_line .= '-';
      }
      $evidence_obj = Bio::PrimarySeq->new(
                      -seq              => $evidence_line,
                      -id               => 0,
  		      -accession_number => $$hit{'hseqname'},
		      -moltype          => $$hit{'moltype'}
		    );
      push @evidence_arr, $evidence_obj;
      $uppercase = 0;	# force next line to start uppercase
    }
    $prev_hseqname = $$hit{'hseqname'};
    $prev_exon = $$hit{'exon'};
  }

  # nucleic acid evidence

  $evidence_obj = Bio::PrimarySeq->new(
                    -seq              => $nucseq,
                    -id               => 0,
     		    -accession_number => $transcript_obj->stable_id,
		    -moltype          => 'dna'
		  );
  push @evidence_arr, $evidence_obj;

  $total_exon_len = 0;
  my @nuc_evidence_arr = ();
  for (my $i = 0; $i <= $#all_exons; $i++) {
    my $start = $all_exons[$i]->start;
    my @features = $all_exons[$i]->each_Supporting_Feature;

    my $last_feat = undef;
    NUC_FEATURE_LOOP:
    foreach my $feature(@features) {
      next NUC_FEATURE_LOOP if ($last_feat
      && ($last_feat->start == $feature->start)
      && ($last_feat->end == $feature->end)
      && ($last_feat->hstart == $feature->hstart)
      && ($last_feat->hseqname eq $feature->hseqname));
      my $fstart = $feature->start;
      my $flen = $feature->end - $fstart + 1;
      my $hlen = $feature->hend -$feature->hstart + 1;
      next NUC_FEATURE_LOOP unless ($flen == $hlen);
      my $hindent;
      my $hit_seq_obj = $self->_get_seq($seqfetcher, \@seqcache,
                                        $feature->hseqname);
      next NUC_FEATURE_LOOP if (! $hit_seq_obj);
      if ($hit_seq_obj->moltype ne "protein") {
        my $hlength = $feature->hend - $feature->hstart + 1;
	my $hseq;
	eval {
          $hseq = substr $hit_seq_obj->seq, $feature->hstart - 1, $hlength;
	};
	next NUC_FEATURE_LOOP if ($@);	# data wrong
        my $strand_wrt_exon = $all_exons[$i]->strand * $feature->strand;
        if ($strand_wrt_exon < 0) {
          my $hseq_obj = Bio::PrimarySeq->new( -seq => $hseq,
                                               -id => 'fake_id',
 		  			       -accession_number => '0',
                                               -moltype => $hit_seq_obj->moltype
                                             );
          $hseq = $hseq_obj->revcom->seq;
        }
        if ($all_exons[$i]->strand > 0) {
          $hindent = $total_exon_len + $fstart - $start;
        } else{
          $hindent = $total_exon_len + $all_exons[$i]->end - $feature->end;
        }
	if ($hindent < 0) {
	  $hindent = 0;
	}
        my %hit_details = ( 'moltype'     => $hit_seq_obj->moltype,
                            'hseqname'    => $feature->hseqname,
                            'hstart'      => $feature->hstart,
                            'hend'        => $feature->hend,
	  		    'hlength'     => $hlength,
			    'hseq'        => $hseq,
			    'hindent'     => $hindent,
			    'hstrand'     => $feature->hstrand,
		            'score'       => $feature->score,
			    'start'       => $start,
			    'end'         => $feature->end,
			    'exon_strand' => $all_exons[$i]->strand,
			    'exon_length' => $all_exons[$i]->end
			                     - $start + 1,
			    'exon'        => $all_exons[$i]
		          );
        push @nuc_evidence_arr, \%hit_details;
      }
      $last_feat = $feature;
    }
    $total_exon_len += $all_exons[$i]->end - $start + 1;
  }

  my @sorted_nuc_evidence_arr = sort {   ( $$a{'hseqname'} =~ /^ENST/
                                           ? -1 : 0 )
			              || ( $$b{'hseqname'} =~ /^ENST/
                                           ? 1 : 0 )
				      || ( $$a{'hseqname'} cmp $$b{'hseqname'} )
				      || ( $$a{'hindent'} <=> $$b{'hindent'} )
 			        } @nuc_evidence_arr;

  $evidence_line = "";
  $uppercase = 1;	# case of sequence for output
  $hit = $sorted_nuc_evidence_arr[0];
  $prev_exon = $$hit{'exon'};
  $prev_hseqname = "-" x 1000;	# fake initial ID
  for (my $i = 0; $i <= $#sorted_nuc_evidence_arr; $i++) {
    my $hit = $sorted_nuc_evidence_arr[$i];
    if ($$hit{'exon'} ne $prev_exon) {
      $uppercase = ! $uppercase;
    }
    my $seq_str = $$hit{'hseq'};
    if ($uppercase) {
      $seq_str = "\U$seq_str";
    } else {
      $seq_str = "\L$seq_str";
    }
    
    if ($$hit{'hseqname'} ne $prev_hseqname) {	# make new evidence line
      $evidence_line = "-" x $nuc_tran_len;
    }

    # splice in the evidence fragment
    substr $evidence_line, $$hit{hindent}, length $seq_str, $seq_str;

    # store if end of evidence line
    if (($i == $#sorted_nuc_evidence_arr) ||
        ($sorted_nuc_evidence_arr[$i+1]{'hseqname'} ne $$hit{'hseqname'}))
    {
      $evidence_line = substr $evidence_line, 0, $nuc_tran_len;
      $evidence_obj = Bio::PrimarySeq->new(
                      -seq              => $evidence_line,
                      -id               => 0,
  		      -accession_number => $$hit{'hseqname'},
		      -moltype          => $$hit{'moltype'}
		    );
      push @evidence_arr, $evidence_obj;
      $uppercase = 0;	# force next line to start uppercase
    }
    $prev_hseqname = $$hit{'hseqname'};
    $prev_exon = $$hit{'exon'};
  }

  return \@evidence_arr;

}

1;

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::SeqIO;
use Bio::EnsEMBL::Gene;

# get a sequence from cache if we have it, otherwise from seqfetcher
sub get_seq {
  die "get_seq: interface fault" if (@_ != 3);
  my ($seqfetcher, $cache_arr_ref, $name) = @_;
  
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

sub get_transcript_nuc {
  die "get_transcript_nuc: interface fault" if (@_ != 1);
  my ($exon_arr_ref) = @_;

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

sub get_transcript_pep {
  die "get_transcript_pep: interface fault" if (@_ != 1);
  my ($exon_pep_arr_ref) = @_;
  
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

# get_aligned_evidence: public subroutine, to be object-ized
# takes a DB adaptor and transcript ID
# returns ref to an array of Bio::PrimarySeq

sub get_aligned_evidence {
  die "get_aligned_evidence: interface fault" if (@_ != 2);
  my ($transcript_id, $db) = @_;

  my $sgp = $db->get_StaticGoldenPathAdaptor;
  my @evidence_arr;	# a reference to this is returned
  my $evidence_obj;

  my $ta = $db->get_TranscriptAdaptor;
  my $ga = $db->get_GeneAdaptor;
  my $ea = $db->get_ExonAdaptor;
  my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher->new();

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

  # protein evidence

  my @utr_free_exons = $transcript_obj->translateable_exons;
  my @exon_peps = ();
  foreach my $exon (@utr_free_exons) {
    if (($exon->end - $exon->start + 1) >= 3) {
      push @exon_peps, $exon->translate()->seq;
    }
  }

  my @all_exons = $transcript_obj->get_all_Exons;

  # record whether each exon is totally UTR, simply a frameshift (hence
  # not to be translated), or at least partly translating; if
  # translating, store the pep sequence string
  my @exon_status = ();
  for (my $i = 0; $i <= $#all_exons; $i++) {
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
    }
  }

  # create an array of all exons that at least partly translate
  my @exons_to_display = ();
  for (my $i = 0; $i <= $#all_exons; $i++) {
    if ($exon_status[$i] eq 'translating') {
      push @exons_to_display, $all_exons[$i];
    }
  }

  $evidence_obj = Bio::PrimarySeq->new(
                    -seq              =>
	              get_transcript_pep(\@exon_peps),
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
    eval {
      print STDERR "exon start end = " .$exons_to_display[$i]->start
      ." ". $exons_to_display[$i]->end . "\n";
    };
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
      my $hit_seq_obj = get_seq($seqfetcher, \@seqcache, $feature->hseqname);
      next PEP_FEATURE_LOOP if (! $hit_seq_obj);
      if ($hit_seq_obj->moltype eq "protein") {
        my $hlength = $feature->hend - $feature->hstart + 1;
        my $length = $feature->end - $feature->start + 1;
        next PEP_FEATURE_LOOP unless ($length == (3 * $hlength));
        printf STDERR "XXX feature length = $length\n";
        my $hseq = substr $hit_seq_obj->seq, $feature->hstart - 1, $hlength;
        my $hindent;
        if ($exons_to_display[$i]->strand > 0) {
          print STDERR "XXX hindent = ($total_exon_len + " . $fstart
 	               . " - " . $start .") / 3";
          $hindent = ($total_exon_len + $fstart - $start) / 3;
  	  print STDERR " = $hindent ( before trunc)\n";
        } else {
          print STDERR "XXX hindent = ($total_exon_len + "
	  . $exons_to_display[$i]->end . " - ". $feature->end . ") /3";
          $hindent = ($total_exon_len + $exons_to_display[$i]->end
	             - $feature->end) / 3;
	  print STDERR " = $hindent (before trunc) \n";
        }
        $hindent = int $hindent;
        if ($i == 0) {
          $hindent -= $utr_region_in_first_exon;
        }
	if ($hindent < 0) {
	  print STDERR "changing hindent $hindent to 0\n";
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
    $total_exon_len += 3 * length($exon_peps[$i]);
  }

  my @sorted_pep_evidence_arr = sort {    $$a{'score'}    <=> $$b{'score'}
                                   || $$a{'hseqname'} cmp $$b{'hseqname'}
 			        } @pep_evidence_arr;

  my $evidence_line = "";
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
      $seq_str = "\L$seq_str";
    }
    print STDERR "XXX ". $$hit{'hseqname'}. " $seq_str\n";
    if ($$hit{'hseqname'} ne $prev_hseqname) {	# not same source as previous
      if ($evidence_line ne "") {
        print STDERR "$evidence_line\n";
	$evidence_obj = Bio::PrimarySeq->new(
	                  -seq              => $evidence_line,
	                  -id               => 0,
			  -accession_number => $$hit{'hseqname'},
			  -moltype          => $$hit{'moltype'}
			);
	push @evidence_arr, $evidence_obj;
      }
      printf STDERR "%-23s", $$hit{'hseqname'}, $$hit{'score'} . "\n";
      $evidence_line = "-" x $$hit{'hindent'};
      print STDERR "XXX evidence_line |$evidence_line|\n";
      print STDERR "XXX adding |$seq_str| to |$evidence_line|\n";
      $evidence_line .= $seq_str;
      print STDERR "XXX evidence_line |$evidence_line|\n";
    } else {	# same source as previous
      while (length($evidence_line) > ($$hit{'hindent'})) {
        chop $evidence_line;
        print STDERR "XXX evidence_line |$evidence_line|\n";
      }
      while (length($evidence_line) < ($$hit{'hindent'})) {
        $evidence_line .= "-";
        print STDERR "XXX evidence_line |$evidence_line|\n";
      }
      print STDERR "XXX adding |$seq_str| to |$evidence_line|\n";
      $evidence_line .= $seq_str;
      print STDERR "XXX evidence_line |$evidence_line|\n";
    }
    $prev_hseqname = $$hit{'hseqname'};
  }
  print STDERR "$evidence_line\n";
  print STDERR "\n";
  if (length($evidence_line) > 0) {
    $evidence_obj = Bio::PrimarySeq->new(
                      -seq              => $evidence_line,
                      -id               => 0,
  		      -accession_number => $$hit{'hseqname'},
		      -moltype          => $$hit{'moltype'}
		    );
    push @evidence_arr, $evidence_obj;
  }

  # nucleic acid evidence

  $evidence_obj = Bio::PrimarySeq->new(
                    -seq              =>
	              get_transcript_nuc(\@all_exons),
                    -id               => 0,
     		    -accession_number => $transcript_obj->stable_id,
		    -moltype          => 'dna'
		  );
  push @evidence_arr, $evidence_obj;

  $total_exon_len = 0;
  my @nuc_evidence_arr = ();
  for (my $i = 0; $i <= $#all_exons; $i++) {
    eval {
      print STDERR "exon start end = " .$all_exons[$i]->start
      ." ". $all_exons[$i]->end . "\n";
    };
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
      my $hit_seq_obj = get_seq($seqfetcher, \@seqcache, $feature->hseqname);
      next NUC_FEATURE_LOOP if (! $hit_seq_obj);
      if ($hit_seq_obj->moltype ne "protein") {
        my $hlength = $feature->hend - $feature->hstart + 1;
        my $hseq = substr $hit_seq_obj->seq, $feature->hstart - 1, $hlength;
        #my $feature_start_overhang = $start - $fstart;
        #my $feature_end_overhang = $feature->end - $all_exons[$i]->end;
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
          #$hindent = $total_exon_len - $feature_start_overhang;
          $hindent = $total_exon_len + $fstart - $start;
        } else{
          $hindent = $total_exon_len + $all_exons[$i]->end - $feature->end;
        }
        #if ($total_exon_len > 0) {
        #  $hindent -= 1;
        #}
	if ($hindent < 0) {
	  print STDERR "changing hindent $hindent to 0\n";
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
    print STDERR "XXX total_exon_len was $total_exon_len\n";
    $total_exon_len += $all_exons[$i]->end - $start + 1;
    print STDERR "XXXXXXX total_exon_len now is " . $total_exon_len . "\n";
  }

  my @sorted_nuc_evidence_arr = sort {    $$a{'score'}    <=> $$b{'score'}
                                       || $$a{'hseqname'} cmp $$b{'hseqname'}
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
    if ($$hit{'hseqname'} ne $prev_hseqname) {	# not same source as previous
      if ($evidence_line ne "") {
        print STDERR "$evidence_line\n";
	$evidence_obj = Bio::PrimarySeq->new(
	                  -seq              => $evidence_line,
	                  -id               => 0,
			  -accession_number => $$hit{'hseqname'},
			  -moltype          => $$hit{'moltype'}
			);
	push @evidence_arr, $evidence_obj;
      }
      printf STDERR "%-23s", $$hit{'hseqname'};
      $evidence_line = "-" x $$hit{'hindent'};
      $evidence_line .= $seq_str;
    } else {	# same source as previous
      while (length($evidence_line) > ($$hit{'hindent'})) {
        chop $evidence_line;
      }
      while (length($evidence_line) < ($$hit{'hindent'})) {
        $evidence_line .= "-";
      }
      $evidence_line .= $seq_str;
    }
    $prev_hseqname = $$hit{'hseqname'};
  }
  print STDERR "$evidence_line\n";
  print STDERR "\n";
  if (length($evidence_line) > 0) {
    $evidence_obj = Bio::PrimarySeq->new(
                      -seq              => $evidence_line,
                      -id               => 0,
    		      -accession_number => $$hit{'hseqname'},
		      -moltype          => $$hit{'moltype'}
		    );
  push @evidence_arr, $evidence_obj;
  }

  return \@evidence_arr;

}

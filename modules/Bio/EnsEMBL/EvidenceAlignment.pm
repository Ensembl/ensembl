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
'-' as required. Proteins are given before DNA. Coordinates in the
database are used to recreate the alignment, and must be correct or some
data may not be displayed.

=head1 CONTACT

Ensembl: ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::EvidenceAlignment;

# modify placement by adding the following to the genomic start/end
use constant MINUS_STRAND_HACK_BP => -1;
use constant PLUS_STRAND_HACK_BP  => +1;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
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
    $obj->{db_adaptor} = $value;
  }
  return $obj->{db_adaptor};
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
    $obj->{transcript_id} = $value;
  }
  return $obj->{transcript_id};
}

=head2 get_features

    Title   :   get_features
    Usage   :   $ea->get_features($transcript_obj, $vc);
    Function:   use SGP adaptor supplied to get evidence off a VC
		of the transcript supplied; return an array
		of featurepairs; features not falling witin exons
		are cut

=cut

sub get_features {
  my ($self, $transcript_obj, $vc) = @_;
  $self->throw('interface fault') if (@_ != 3);
    
  my @exons = $transcript_obj->get_all_Exons;
  my $strand = $exons[0]->strand;
  print STDERR "transcript ", $transcript_obj->stable_id, "strand ", $strand, "\n";
  my @all_features = $vc->get_all_SimilarityFeatures;
  my @features = ();
  FEATURE_LOOP:
  foreach my $feature (@all_features) {
    if ($feature->strand == $strand) {
      foreach my $exon (@exons) {
        if ($feature->start >= $exon->start and $feature->end <= $exon->end) {
	  push @features, $feature;
	  next FEATURE_LOOP;
	}
      }
    }
  }
  return @features;
}

=head2 pad_pep_str

    Title   :   pad_pep_str
    Usage   :   $padded_pep_seq = $ea->pad_pep_str($pep_str);
    Function:   make a protein sequence the same length as the
                corresponding DNA sequence by adding '--' after
		each amino acid
    Returns :   string

=cut

sub pad_pep_str {
  my ($self, $original) = @_;
  $self->throw('interface fault') if (@_ != 2);

  my @amino_acids = split //, $original;
  my $padded = '';
  foreach my $amino_acid (@amino_acids) {
    $padded .= $amino_acid . '--';
  }
  return $padded;
}

=head2 fetch_alignment

    Title   :   fetch_alignment
    Usage   :   my $seq_arr_ref = $ea->fetch_alignment;
    Function:   gets aligned transcript and evidence
    Returns :   reference to array of Bio::PrimarySeq

=cut

sub fetch_alignment {
  my ($self) = @_;
  $self->throw('interface fault') if (@_ != 1);
  $self->throw('must have a transcript stable ID and a DB adaptor object')
    unless ($self->transcriptid && $self->dbadaptor);
  return $self->_get_aligned_evidence($self->transcriptid, $self->dbadaptor);
}

# sort an evidence array
sub _evidence_sort {
  my ($self, $evidence_arr_ref) = @_;
  $self->throw('interface fault') if (@_ != 2);

  # ENST... first. Then alphabetically by name, followed by indent
  my @sorted_evidence_arr = sort {    ( $$a{hseqname} =~ /^ENST/ ? -1 : 0 )
			           || ( $$b{hseqname} =~ /^ENST/ ?  1 : 0 )
				   || ( $$a{hseqname} cmp $$b{hseqname} )
				   || ( $$a{hindent}  <=> $$b{hindent} )
 			         } @$evidence_arr_ref;
  return \@sorted_evidence_arr;
}

# get cDNA
sub _get_transcript_nuc {
  my ($self, $exon_arr_ref) = @_;
  $self->throw('interface fault') if (@_ != 2);

  my $retval = '';
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

# get all hit sequences by aggregated use of pfetch
sub _get_hits {
  my ($self, $features_arr_ref) = @_;
  $self->throw('interface fault') if (@_ != 2);

  my $pfetch = '/usr/local/pubseq/bin/pfetch';	# executable
  my $clump_size = 1000;	# number of sequences to fetch at once
  my %hits_hash = ();
  my @hseqnames = ();
  for (my $i = 0; $i < @$features_arr_ref; $i++) {
    push @hseqnames, $$features_arr_ref[$i]->hseqname;
    if (($i % $clump_size) == 0) {
      open (PFETCH_IN, "$pfetch -q @hseqnames |")
        or $self->throw("error running pfetch");
      my $seq_no = 0;
      while (<PFETCH_IN>) {
        chomp;
	my $seq_obj;
	if ($_ ne "no match") {
	   my $seq_obj = Bio::Seq->new( -seq => $_,
	                                -id  => 'fake_id',
				        -accession_number =>$hseqnames[$seq_no]
				      );
 	  $hits_hash{$hseqnames[$seq_no]} = $seq_obj;
	}
	$seq_no++;
      }
      @hseqnames = ();
    }
  }

  # fetch the non-clump-sized remainder
  open (PFETCH_IN, "$pfetch -q @hseqnames |")
    or $self->throw("error running pfetch");
  my $seq_no = 0;
  while (<PFETCH_IN>) {
    chomp;
    my $seq_obj;
    if ($_ ne "no match") {
      my $seq_obj = Bio::Seq->new( -seq => $_,
                                   -id  => 'fake_id',
				   -accession_number =>$hseqnames[$seq_no]
				 );
      $hits_hash{$hseqnames[$seq_no]} = $seq_obj;
    }
    $seq_no++;
  }

  return \%hits_hash;
}

# _get_aligned_evidence: takes a transcript ID and a DB adaptor
# returns ref to an array of Bio::PrimarySeq

sub _get_aligned_evidence {
  my ($self, $transcript_id, $db) = @_;
  $self->throw('interface fault') if (@_ != 3);

  my $sgp = $db->get_StaticGoldenPathAdaptor;
  my @evidence_arr;	# a reference to this is returned
  my $evidence_obj;

  my $ta = $db->get_TranscriptAdaptor;
  my $ga = $db->get_GeneAdaptor;
  my $ea = $db->get_ExonAdaptor;

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

  my @all_exons = $transcript_obj->get_all_Exons;
  my @features = $self->get_features($transcript_obj, $vc);
  my $hits_hash_ref = $self->_get_hits(\@features);
  my $translation = $transcript_obj->translate->seq;
  my $nucseq = $self->_get_transcript_nuc(\@all_exons);
  my $cdna_len_bp = length $nucseq;

  # protein evidence

  my @utr_free_exons = $transcript_obj->translateable_exons;

  # record whether each exon is totally UTR, simply a frameshift (hence
  # not to be translated), or at least partly translating
  my @exon_status = ();
  for (my $i = 0; $i < @all_exons; $i++) {
    $exon_status[$i] = 'utr';
  }
  for (my $i = 0; $i < @all_exons; $i++) {
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
  my $total_5prime_utr_len = 0;
  my $fiveprime_utr_in_first_exon = 0;
  UTR_OUTER:
  for (my $i = 0; $i < @all_exons; $i++) {
    my $exon_i = $all_exons[$i];
    if ($exon_status[$i] ne 'translating') {
      $total_5prime_utr_len += $exon_i->end - $exon_i->start + 1;
    } else {	# 1st translating exon, could contain part of 5' UTR
      foreach my $utr_free_exon (@utr_free_exons) {
        if (($utr_free_exon->start >= $exon_i->start) 
        && ($utr_free_exon->end <= $exon_i->end))
	{
	  if ($exon_i->strand > 0) {
            $fiveprime_utr_in_first_exon =  $utr_free_exon->start
	                                  - $exon_i->start;
          } else {
	    $fiveprime_utr_in_first_exon = $exon_i->end - $utr_free_exon->end;
          }
	  last UTR_OUTER;
	}
      }
    }
  }
  if ($fiveprime_utr_in_first_exon < 0) {
    $fiveprime_utr_in_first_exon = 0;	# disaster recovery
  }
  $total_5prime_utr_len += $fiveprime_utr_in_first_exon;

  # translation itself forms our first row of 'evidence'
  my $evidence_line = $translation;
  $evidence_line = $self->pad_pep_str($evidence_line);
  # adjust for 5' and 3' UTRs
  $evidence_line = ('-' x $total_5prime_utr_len) . $evidence_line;
  while (length($evidence_line) < $cdna_len_bp) {
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
  my $total_exon_len = 0;
  my @pep_evidence_arr = ();
  my $last_feat = undef;
  foreach my $exon (@all_exons) {
    PEP_FEATURE_LOOP:
    foreach my $feature(@features) {
      next PEP_FEATURE_LOOP	# unless feature falls within this exon
        unless (($feature->start >= $exon->start)
	     && ($feature->end <= $exon->end));
      next PEP_FEATURE_LOOP	# if feature is a duplicate of the last
        if ($last_feat
        && ($last_feat->start  == $feature->start)
        && ($last_feat->end    == $feature->end)
        && ($last_feat->hstart == $feature->hstart)
        && ($last_feat->hseqname eq $feature->hseqname));
      my $hit_seq_obj = $$hits_hash_ref{$feature->hseqname};
      next PEP_FEATURE_LOOP 
        unless ($hit_seq_obj && ($hit_seq_obj->moltype eq 'protein'));
      my $hlen = $feature->hend - $feature->hstart + 1;
      my $flen = $feature->end - $feature->start + 1;
      next PEP_FEATURE_LOOP unless ($flen == 3 * $hlen);
      if (($feature->hstart - 1 < 0) || ($feature->hstart - 1 + $hlen
	  > length $hit_seq_obj->seq))
      {
        print STDERR 'XXX pep hit coords out of range for ',
          $feature->hseqname, ': hstart ', $feature->hstart, ', hend ',
          $feature->hend, ', max length ', length $hit_seq_obj->seq, "\n";
        next PEP_FEATURE_LOOP;
      }
      my $hseq = substr $hit_seq_obj->seq, $feature->hstart - 1, $hlen;
      $hseq = $self->pad_pep_str($hseq);
      my $hindent_bp;
      if ($exon->strand > 0) {
        $hindent_bp =   $total_exon_len + $feature->start - $exon->start
	              + PLUS_STRAND_HACK_BP;
      } else {
        $hindent_bp =   $total_exon_len + $exon->end - $feature->end
	              + MINUS_STRAND_HACK_BP;
      }
      if ($hindent_bp < 0) {
        $hindent_bp = 0;	# disaster recovery
      }
      my %hit_details = ( 'moltype'  => $hit_seq_obj->moltype,
                          'hseqname' => $feature->hseqname,
			  'hseq'     => $hseq,
			  'hindent'  => $hindent_bp,
                          'exon'     => $exon
			);
      push @pep_evidence_arr, \%hit_details;
      $last_feat = $feature;
    }
    $total_exon_len += $exon->end - $exon->start + 1;
  }

  my @sorted_pep_evidence_arr = @{$self->_evidence_sort(\@pep_evidence_arr)};

  $evidence_line = '';
  my $prev_hseqname = '?' x 1000;	# fake initial ID
  for (my $i = 0; $i < @sorted_pep_evidence_arr; $i++) {
    my $hit = $sorted_pep_evidence_arr[$i];
    if ($$hit{hseqname} ne $prev_hseqname) {	# make new evidence line
      $evidence_line = '-' x $cdna_len_bp;
    }

    # splice in the evidence fragment
    my $hseqlen = length $$hit{hseq};
    next if (($$hit{hindent} < $total_5prime_utr_len)
          || ($$hit{hindent} + $hseqlen > $cdna_len_bp));

    substr $evidence_line, $$hit{hindent}, $hseqlen, $$hit{hseq};

    # store if end of evidence line
    if (($i == $#sorted_pep_evidence_arr)
     || ($sorted_pep_evidence_arr[$i+1]{hseqname} ne $$hit{hseqname}))
    {
      $evidence_obj = Bio::PrimarySeq->new(
                      -seq              => $evidence_line,
                      -id               => 0,
  		      -accession_number => $$hit{hseqname},
		      -moltype          => $$hit{moltype}
		    );
      push @evidence_arr, $evidence_obj;
    }
    $prev_hseqname = $$hit{hseqname};
  }

  # nucleic acid evidence

  # cDNA itself forms first row of nucleic acid evidence
  $evidence_obj = Bio::PrimarySeq->new(
                    -seq              => $nucseq,
                    -id               => 0,
     		    -accession_number => $transcript_obj->stable_id,
		    -moltype          => 'dna'
		  );
  push @evidence_arr, $evidence_obj;

  $total_exon_len = 0;
  my @nuc_evidence_arr = ();
  $last_feat = undef;
  foreach my $exon (@all_exons) {
    NUC_FEATURE_LOOP:
    foreach my $feature(@features) {
      next NUC_FEATURE_LOOP	# unless feature falls within this exon
        unless (($feature->start >= $exon->start)
	     && ($feature->end <= $exon->end));
      next NUC_FEATURE_LOOP	# if feature is a duplicate of the last
        if ($last_feat
        && ($last_feat->start  == $feature->start)
        && ($last_feat->end    == $feature->end)
        && ($last_feat->hstart == $feature->hstart)
        && ($last_feat->hseqname eq $feature->hseqname));
      my $hit_seq_obj = $$hits_hash_ref{$feature->hseqname};
      next NUC_FEATURE_LOOP 
        unless ($hit_seq_obj && ($hit_seq_obj->moltype ne 'protein'));
      my $hlen = $feature->hend - $feature->hstart + 1;
      my $flen = $feature->end - $feature->start + 1;
      next NUC_FEATURE_LOOP unless ($flen == $hlen);
      if (($feature->hstart - 1 < 0) || ($feature->hstart - 1 + $hlen
	  > length $hit_seq_obj->seq))
      {
        print STDERR 'XXX nuc hit coords out of range for ',
          $feature->hseqname, ': hstart ', $feature->hstart, ', hend ',
          $feature->hend, ', max length ', length $hit_seq_obj->seq, "\n";
        next NUC_FEATURE_LOOP;
      }
      my $hseq = substr $hit_seq_obj->seq, $feature->hstart - 1, $hlen;
      my $strand_wrt_exon = $exon->strand * $feature->strand;
      if ($strand_wrt_exon < 0) {	# reverse-compliment the hit
        my $hseq_obj_tmp = Bio::PrimarySeq->new(
	                                    -seq => $hseq,
                                            -id => 'fake_id',
                                            -accession_number => '0',
                                            -moltype => $hit_seq_obj->moltype
                                           );
        $hseq = $hseq_obj_tmp->revcom->seq;
      }
      my $hindent_bp;
      if ($exon->strand > 0) {
        $hindent_bp =   $total_exon_len + $feature->start - $exon->start
	              + PLUS_STRAND_HACK_BP;
      } else{
        $hindent_bp =   $total_exon_len + $exon->end - $feature->end
	              + MINUS_STRAND_HACK_BP;
      }
      if ($hindent_bp < 0) {
        $hindent_bp = 0;	# disaster recovery
      }
      my %hit_details = ( 'moltype'     => $hit_seq_obj->moltype,
                          'hseqname'    => $feature->hseqname,
			  'hseq'        => $hseq,
			  'hindent'     => $hindent_bp,
			  'exon'        => $exon
		          );
      push @nuc_evidence_arr, \%hit_details;
      $last_feat = $feature;
    }
    $total_exon_len += $exon->end - $exon->start + 1;
  }

  my @sorted_nuc_evidence_arr = @{$self->_evidence_sort(\@nuc_evidence_arr)};

  $evidence_line = '';
  my $uppercase = 1;	# case of sequence for output
  my $hit = $sorted_nuc_evidence_arr[0];
  my $prev_exon = $$hit{exon};
  $prev_hseqname = '-' x 1000;	# fake initial ID
  for (my $i = 0; $i < @sorted_nuc_evidence_arr; $i++) {
    my $hit = $sorted_nuc_evidence_arr[$i];
    if ($$hit{exon} ne $prev_exon) {
      $uppercase = ! $uppercase;
    }
    my $hseq_str = $$hit{hseq};
    if ($uppercase) {
      $hseq_str = "\U$hseq_str";
    } else {
      $hseq_str = "\L$hseq_str";
    }
    if ($$hit{hseqname} ne $prev_hseqname) {	# make new evidence line
      $evidence_line = '-' x $cdna_len_bp;
    }

    # splice in the evidence fragment
    my $hseqlen = length $$hit{hseq};
    next unless ($$hit{hindent} + $hseqlen <= $cdna_len_bp);
    substr $evidence_line, $$hit{hindent}, $hseqlen, $hseq_str;

    # store if end of evidence line
    if (($i == $#sorted_nuc_evidence_arr)
     || ($sorted_nuc_evidence_arr[$i+1]{hseqname} ne $$hit{hseqname}))
    {
      $evidence_obj = Bio::PrimarySeq->new(
                      -seq              => $evidence_line,
                      -id               => 0,
  		      -accession_number => $$hit{hseqname},
		      -moltype          => $$hit{moltype}
		    );
      push @evidence_arr, $evidence_obj;
      $uppercase = 0;	# so next line starts uppercase
    }
    $prev_hseqname = $$hit{hseqname};
    $prev_exon = $$hit{exon};
  }

  # remove blank evidence lines
  my @filtered_evidence_arr = ();
  foreach my $evidence_line (@evidence_arr) {
    push @filtered_evidence_arr, $evidence_line
      if ($$evidence_line{seq} =~ /[^-]/);
  }
  
  return \@filtered_evidence_arr;

}

1;

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
 my $alignment_arr_ref = $ea->fetch_alignment;
 $ea->transcriptid($other_tr_stable_id);
 my $other_alignment_arr_ref = $ea->fetch_alignment;
 $ea->contigid($contig_stable_id);
 my $contig_based_alignment_arr_ref = $ea->fetch_alignment;

=head1 DESCRIPTION

Gives an alignment of either a transcript and its evidence or a raw
contig and its similarity features. Coordinates in the database are
used to recreate the alignment, and must be correct or some data may
not be displayed or may be displayed incorrectly.

=head1 CONTACT

Ensembl: ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::EvidenceAlignment;

# pfetch binary - just pick it up from PATH
use constant PFETCH => 'pfetch';

# modify placement on VC by adding the following to genomic start/end
use constant VC_MINUS_STRAND_HACK_BP => -1;
use constant VC_PLUS_STRAND_HACK_BP  => +1;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Gene;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

=head2 new

    Title   :   new
    Usage   :   my $tr_ea   = Bio::EnsEMBL::EvidenceAlignment->new(
                                -DBADAPTOR    => $dba,
                                -TRANSCRIPTID => $transcript_stable_id);
		my $cont_ea = Bio::EnsEMBL::EvidenceAlignment->new(
                                -DBADAPTOR => $dba,
                                -CONTIGID  => $contig_stable_id);
    Function:   Initialises EvidenceAlignment object
    Returns :   An EvidenceAlignment object
    Args    :   Database adaptor object and an ID string (-CONTIGID
                with contig ID or -TRANSCRIPTID with transcript
		stable ID).

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($transcriptid, $contigid, $dbadaptor) = $self->_rearrange(
                                                ['TRANSCRIPTID',
						 'CONTIGID',
                                                 'DBADAPTOR'],
						@args);
  if (defined $transcriptid and defined $contigid) {
    $self->throw("may have a transcript ID or a contig ID but not both");
  }
  if ($transcriptid) {
    $self->transcriptid($transcriptid);
  }
  if ($contigid) {
    $self->contigid($contigid);
  }
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
    $obj->{evidencealignment_db_adaptor} = $value;
  }
  return $obj->{evidencealignment_db_adaptor};
}

=head2 transcriptid

    Title   :   transcriptid
    Usage   :   $ea->transcriptid($transcript_stable_id);
    Function:   get/set for transcript stable id string (setting
                this also unsets contigid)

=cut

sub transcriptid {
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{evidencealignment_transcript_id} = $value;
    undef $obj->{evidencealignment_contig_id};
  }
  return $obj->{evidencealignment_transcript_id};
}

=head2 contigid

    Title   :   contigid
    Usage   :   $ea->contigid($contig_stable_id);
    Function:   get/set for contig stable id string (setting this
                also unsets transcriptid)

=cut

sub contigid {
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{evidencealignment_contig_id} = $value;
    undef $obj->{evidencealignment_transcript_id};
  }
  return $obj->{evidencealignment_contig_id};
}

=head2 _get_features_from_transcript

    Title   :   _get_features_from_transcript
    Usage   :   $ea->_get_features_from_transcript($transcript_obj, $vc);
    Function:   use SGP adaptor supplied to get evidence off a VC
		of the transcript supplied; features not falling witin
		exons are cut
    Returns :   array of featurepairs

=cut

sub _get_features_from_transcript {
  my ($self, $transcript_obj, $vc) = @_;
  $self->throw('interface fault') if (@_ != 3);
    
  my @exons = $transcript_obj->get_all_Exons;
  my $strand = $exons[0]->strand;
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

=head2 _get_features_from_rawcontig

    Title   :   _get_features_from_rawcontig
    Usage   :   $ea->_get_features_from_rawcontig($rc_obj, $strand);
    Function:   Get features off specified strand of the raw contig
                supplied
    Returns :   array of featurepairs

=cut

sub _get_features_from_rawcontig {
  my ($self, $rawcontig_obj, $strand) = @_;
  $self->throw('interface fault') if (@_ != 3);

  my @all_features = $rawcontig_obj->get_all_SimilarityFeatures;
  my @features = ();
  foreach my $feature (@all_features) {
    if ($feature->strand == $strand) {
      my $tmp;
      eval {
        $tmp = $feature->hseqname;
      };
      if ($@) {
        $self->warn("feature $feature has no hseqname method");
      } else {
        push @features, $feature;
      }
    }
  }
  return @features;
}

=head2 _pad_pep_str

    Title   :   _pad_pep_str
    Usage   :   $padded_pep_seq = $ea->_pad_pep_str($pep_str);
    Function:   make a protein sequence the same length as the
                corresponding DNA sequence by adding '--' after
		each amino acid
    Returns :   string

=cut

sub _pad_pep_str {
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
    Function:   gets transcript or raw contig and corresponding
                evidence or similarity features; for raw
		contigs, these are displayed for the forward
		strand followed by the reverse strand
    Returns :   reference to array of Bio::PrimarySeq

=cut

sub fetch_alignment {
  my ($self) = @_;
  $self->throw('interface fault') if (@_ != 1);
  $self->throw('must have a stable ID and a DB adaptor object')
    unless (($self->transcriptid || $self->contigid) && $self->dbadaptor);
  if ($self->transcriptid) {
    return $self->_get_aligned_evidence_for_transcript($self->transcriptid,
                                                       $self->dbadaptor);
  } elsif ($self->contigid) {
    my $plus_strand_alignment  = $self->_get_aligned_features_for_contig(
                                 $self->contigid, $self->dbadaptor, 1);
    my $minus_strand_alignment = $self->_get_aligned_features_for_contig(
                                 $self->contigid, $self->dbadaptor, -1);
    my $all_alignments = $plus_strand_alignment;
    foreach my $line (@$minus_strand_alignment) {
      push @$all_alignments, $line;
    }
    return $all_alignments;
  }
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
  foreach my $exon (@$exon_arr_ref) {
    $retval .= $exon->seq->seq;
  }
  return $retval;
}

# get all hit sequences by aggregated use of pfetch
sub _get_hits {
  my ($self, $features_arr_ref) = @_;
  $self->throw('interface fault') if (@_ != 2);

  my $pfetch = PFETCH;		# executable
  my $clump_size = 1000;	# number of sequences to fetch at once
  my %hits_hash = ();
  my @hseqnames = ();
  for (my $i = 0; $i < @$features_arr_ref; $i++) {
    my $hseqname = $$features_arr_ref[$i]->hseqname;
    if (! exists $hits_hash{$hseqname})
    {
      push @hseqnames, $$features_arr_ref[$i]->hseqname;
      if ((@hseqnames % $clump_size) == 0) {
        open (EVIDENCEALIGNMENT_PFETCH_IN_FH, "$pfetch -q @hseqnames |")
          or $self->throw("error running pfetch");
        my $seq_no = 0;
        while (<EVIDENCEALIGNMENT_PFETCH_IN_FH>) {
          chomp;
	  my $seq_obj;
	  if ($_ ne "no match") {
	     my $seq_obj = Bio::Seq->new(
	                     -seq => $_,
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
  }

  # fetch the non-clump-sized remainder
  open (EVIDENCEALIGNMENT_PFETCH_IN_FH, "$pfetch -q @hseqnames |")
    or $self->throw("error running pfetch");
  my $seq_no = 0;
  while (<EVIDENCEALIGNMENT_PFETCH_IN_FH>) {
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

# _get_aligned_features_for_contig: takes a contig ID, a DB adaptor
# and strand
# returns ref to an array of Bio::PrimarySeq

sub _get_aligned_features_for_contig {
  my ($self, $contig_id, $db, $strand) = @_;
  $self->throw('interface fault') if (@_ != 4);

  my @evidence_arr;	# a reference to this is returned
  my $evidence_obj;

  # get contig
  my $contig_obj = $db->get_Contig($contig_id);

  my @features = $self->_get_features_from_rawcontig($contig_obj, $strand);
  my $hits_hash_ref = $self->_get_hits(\@features);
  my $nucseq_obj = $contig_obj->primary_seq;
  if ($strand < 0) {
    $nucseq_obj = $nucseq_obj->revcom;
  }
  my $nucseq_str = $nucseq_obj->seq;
  my $dna_len_bp = length $nucseq_str;
  my @translations = ();
  for (my $i = 0; $i < 3; $i++) {
    push @translations, $nucseq_obj->translate(undef, undef, $i);
  }

  # protein evidence

  # translations themselves form our first three rows of 'evidence'
  for (my $i = 0; $i < 3; $i++) {
    my $evidence_line = $translations[$i]->seq;
    $evidence_line = $self->_pad_pep_str($evidence_line);
    if ($i == 1) {
      $evidence_line = '-' . $evidence_line;
    } elsif ($i == 2) {
      $evidence_line = '--' . $evidence_line;
    }
    while (length $evidence_line < $dna_len_bp) {
      $evidence_line .= '-';
    }
    while (length $evidence_line > $dna_len_bp) {
      chop $evidence_line;
    }
    $evidence_obj = Bio::PrimarySeq->new(
                      -seq              => $evidence_line,
                      -id               => 0,
                      -accession_number => $contig_obj->id,
                      -moltype          => 'protein'
		    );
    push @evidence_arr, $evidence_obj;
  }

  my @seqcache = ();
  my @pep_evidence_arr = ();
  my $last_feat = undef;

  PEP_FEATURE_LOOP:
  foreach my $feature (@features) {
    next PEP_FEATURE_LOOP	# if feature is a duplicate of the last
      if ($last_feat
      && ($last_feat->start  == $feature->start)
      && ($last_feat->end    == $feature->end)
      && ($last_feat->hstart == $feature->hstart)
      && ($last_feat->hseqname eq $feature->hseqname));
    if (($feature->start < 1) || ($feature->end > $dna_len_bp)) {
      $self->warn("genomic coordinates out of range: start " .
        $feature->start . ", end " . $feature->end);
      next PEP_FEATURE_LOOP;
    }
    my $hit_seq_obj = $$hits_hash_ref{$feature->hseqname};
    if (! $hit_seq_obj) {
      $self->warn("couldn't fetch hit sequence " . $feature->hseqname . "\n");
      next PEP_FEATURE_LOOP;
    }
    next PEP_FEATURE_LOOP	# not an error, DNA and protein are mixed
      unless ($hit_seq_obj->moltype eq 'protein');
    my $hlen = $feature->hend - $feature->hstart + 1;
    my $flen = $feature->end - $feature->start + 1;
    if ($flen != 3 * $hlen) {
      $self->warn("genomic length $flen but protein hit length $hlen for hit "
        . $feature->hseqname . "\n");
      next PEP_FEATURE_LOOP;
    }
    if (($feature->hstart - 1 < 0) || ($feature->hstart - 1 + $hlen
      > length $hit_seq_obj->seq))
    {
      $self->warn("hit coordinates out of range: hit " . $feature->hseqname .
        ", hit start " . $feature->hstart . ", hit length $hlen\n");
      next PEP_FEATURE_LOOP;
    }
    my $hseq = substr $hit_seq_obj->seq, $feature->hstart - 1, $hlen;
    $hseq = $self->_pad_pep_str($hseq);
    my $hindent_bp;
    if ($strand > 0) {
      $hindent_bp = $feature->start - 1;
    } else {
      $hindent_bp = $dna_len_bp - $feature->end;
    }
    if ($hindent_bp < 0) {
      $hindent_bp = 0;	# disaster recovery
    }
    my %hit_details = ( 'moltype'  => $hit_seq_obj->moltype,
                        'hseqname' => $feature->hseqname,
                        'hseq'     => $hseq,
                        'hindent'  => $hindent_bp,
                      );
    push @pep_evidence_arr, \%hit_details;
    $last_feat = $feature;
  }

  my @sorted_pep_evidence_arr = @{$self->_evidence_sort(\@pep_evidence_arr)};

  my $evidence_line = '';
  my $prev_hseqname = '-' x 1000;	# fake initial ID
  for (my $i = 0; $i < @sorted_pep_evidence_arr; $i++) {
    my $hit = $sorted_pep_evidence_arr[$i];
    if ($$hit{hseqname} ne $prev_hseqname) {	# make new evidence line
      $evidence_line = '-' x $dna_len_bp;
    }

    # splice in the evidence fragment
    my $hseqlen = length $$hit{hseq};
    next
      if (($$hit{hindent} < 0) || ($$hit{hindent} + $hseqlen > $dna_len_bp));
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

  # DNA itself forms first row of nucleic acid evidence
  $evidence_obj = Bio::PrimarySeq->new(
                    -seq              => $nucseq_str,
                    -id               => 0,
     		    -accession_number => $contig_obj->id,
		    -moltype          => 'dna'
		  );
  push @evidence_arr, $evidence_obj;

  my @nuc_evidence_arr = ();
  $last_feat = undef;
  NUC_FEATURE_LOOP:
  foreach my $feature(@features) {
    next NUC_FEATURE_LOOP	# if feature is a duplicate of the last
      if ($last_feat
      && ($last_feat->start  == $feature->start)
      && ($last_feat->end    == $feature->end)
      && ($last_feat->hstart == $feature->hstart)
      && ($last_feat->hseqname eq $feature->hseqname));
    if (($feature->start < 1) || ($feature->end > $dna_len_bp)) {
      $self->warn("genomic coordinates out of range: start " .
        $feature->start . ", end " . $feature->end);
      next NUC_FEATURE_LOOP;
    }
    my $hit_seq_obj = $$hits_hash_ref{$feature->hseqname};
    if (! $hit_seq_obj) {
      $self->warn("couldn't fetch hit sequence " . $feature->hseqname . "\n");
      next NUC_FEATURE_LOOP;
    }
    next NUC_FEATURE_LOOP	# not an error, DNA and protein are mixed
      unless ($hit_seq_obj->moltype ne 'protein');
    my $hlen = $feature->hend - $feature->hstart + 1;
    my $flen = $feature->end - $feature->start + 1;
    if ($hlen != $flen) {
      $self->warn("genomic length $flen but DNA hit length $hlen for hit "
        . $feature->hseqname . "\n");
      next NUC_FEATURE_LOOP;
    }
    if (($feature->hstart - 1 < 0) || ($feature->hstart - 1 + $hlen
         > length $hit_seq_obj->seq))
    {
      $self->warn("hit coordinates out of range: hit " . $feature->hseqname .
        ", hit start " . $feature->hstart . ", hit length $hlen\n");
      next NUC_FEATURE_LOOP;
    }
    my $hseq = substr $hit_seq_obj->seq, $feature->hstart - 1, $hlen;
    my $strand_wrt_dna = $strand * $feature->strand;
    if ($strand_wrt_dna < 0) {	# reverse-compliment the hit
      my $hseq_obj_tmp = Bio::PrimarySeq->new(
                                          -seq => $hseq,
                                          -id => 'fake_id',
                                          -accession_number => '0',
                                          -moltype => $hit_seq_obj->moltype
                                         );
      $hseq = $hseq_obj_tmp->revcom->seq;
    }
    my $hindent_bp;
    if ($strand > 0) {
      $hindent_bp = $feature->start - 1;
    } else{
      $hindent_bp = $dna_len_bp - $feature->end;
    }
    if ($hindent_bp < 0) {
      $hindent_bp = 0;	# disaster recovery
    }
    my %hit_details = ( 'moltype'     => $hit_seq_obj->moltype,
                        'hseqname'    => $feature->hseqname,
                        'hseq'        => $hseq,
		        'hindent'     => $hindent_bp,
		      );
    push @nuc_evidence_arr, \%hit_details;
    $last_feat = $feature;
  }

  my @sorted_nuc_evidence_arr = @{$self->_evidence_sort(\@nuc_evidence_arr)};

  $evidence_line = '';
  my $hit = $sorted_nuc_evidence_arr[0];
  $prev_hseqname = '-' x 1000;	# fake initial ID
  for (my $i = 0; $i < @sorted_nuc_evidence_arr; $i++) {
    my $hit = $sorted_nuc_evidence_arr[$i];
    my $hseq_str = $$hit{hseq};
    if ($$hit{hseqname} ne $prev_hseqname) {	# make new evidence line
      $evidence_line = '-' x $dna_len_bp;
    }

    # splice in the evidence fragment
    my $hseqlen = length $$hit{hseq};
    next unless ($$hit{hindent} + $hseqlen <= $dna_len_bp);
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
    }
    $prev_hseqname = $$hit{hseqname};
  }

  # remove blank evidence lines

  my @filtered_evidence_arr = ();
  foreach my $evidence_line (@evidence_arr) {
    push @filtered_evidence_arr, $evidence_line
      if ($$evidence_line{seq} =~ /[^-]/);
  }

  # use uppercase
  foreach my $evidence_line (@filtered_evidence_arr) {
    my $seq_str = uc $evidence_line->seq;
    $evidence_line->seq($seq_str);
  }

  return \@filtered_evidence_arr;
}

# _get_aligned_evidence_for_transcript: takes a transcript ID and
# a DB adaptor
# returns ref to an array of Bio::PrimarySeq

sub _get_aligned_evidence_for_transcript {
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
  my @features = $self->_get_features_from_transcript($transcript_obj, $vc);
  my $hits_hash_ref = $self->_get_hits(\@features);
  my $translation = $transcript_obj->translate->seq;
  my $nucseq_str = $self->_get_transcript_nuc(\@all_exons);
  my $cdna_len_bp = length $nucseq_str;

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
  $evidence_line = $self->_pad_pep_str($evidence_line);
  # get 3' UTR length and adjust evidence line for 5' and 3' UTRs
  $evidence_line = ('-' x $total_5prime_utr_len) . $evidence_line;
  my $total_3prime_utr_len = $cdna_len_bp - length($evidence_line);
  if ($total_3prime_utr_len < 0) {
    $total_3prime_utr_len = 0;	# disaster recovery
  }
  $evidence_line .= '-' x $total_3prime_utr_len;
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
      if (! $hit_seq_obj) {
        $self->warn("couldn't fetch hit sequence " . $feature->hseqname ."\n");
        next PEP_FEATURE_LOOP;
      }
      next PEP_FEATURE_LOOP	# not an error, DNA and protein are mixed
        unless ($hit_seq_obj->moltype eq 'protein');
      my $hlen = $feature->hend - $feature->hstart + 1;
      my $flen = $feature->end - $feature->start + 1;
      if ($flen != 3 * $hlen) {
        $self->warn("genomic length $flen but protein hit length $hlen for hit "
        . $feature->hseqname . "\n");
        next PEP_FEATURE_LOOP;
      }
      if (($feature->hstart - 1 < 0) || ($feature->hstart - 1 + $hlen
	  > length $hit_seq_obj->seq))
      {
        $self->warn("hit coordinates out of range: hit " . $feature->hseqname .
          ", hit start " . $feature->hstart . ", hit length $hlen\n");
        next PEP_FEATURE_LOOP;
      }
      my $hseq = substr $hit_seq_obj->seq, $feature->hstart - 1, $hlen;
      $hseq = $self->_pad_pep_str($hseq);
      my $hindent_bp;
      if ($exon->strand > 0) {
        $hindent_bp =   $total_exon_len + $feature->start - $exon->start
	              + VC_PLUS_STRAND_HACK_BP;
      } else {
        $hindent_bp =   $total_exon_len + $exon->end - $feature->end
	              + VC_MINUS_STRAND_HACK_BP;
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
  my $prev_hseqname = '-' x 1000;	# fake initial ID
  for (my $i = 0; $i < @sorted_pep_evidence_arr; $i++) {
    my $hit = $sorted_pep_evidence_arr[$i];
    if ($$hit{hseqname} ne $prev_hseqname) {	# make new evidence line
      $evidence_line = '-' x $cdna_len_bp;
    }

    # splice in the evidence fragment
    my $hseqlen = length $$hit{hseq};
    next if (($$hit{hindent} < $total_5prime_utr_len)
          || ($$hit{hindent} + $hseqlen
	     > ($cdna_len_bp - $total_3prime_utr_len)));
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
                    -seq              => $nucseq_str,
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
      if (! $hit_seq_obj) {
        $self->warn("couldn't fetch hit sequence " . $feature->hseqname ."\n");
      }
      next NUC_FEATURE_LOOP	# not an error, DNA and protein are mixed
        unless ($hit_seq_obj->moltype ne 'protein');
      my $hlen = $feature->hend - $feature->hstart + 1;
      my $flen = $feature->end - $feature->start + 1;
      if ($hlen != $flen) {
        $self->warn("genomic length $flen but DNA hit length $hlen for hit "
          . $feature->hseqname . "\n");
        next NUC_FEATURE_LOOP;
      }
      if (($feature->hstart - 1 < 0) || ($feature->hstart - 1 + $hlen
	  > length $hit_seq_obj->seq))
      {
        $self->warn("hit coordinates out of range: hit " . $feature->hseqname .
          ", hit start " . $feature->hstart . ", hit length $hlen\n");
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
	              + VC_PLUS_STRAND_HACK_BP;
      } else{
        $hindent_bp =   $total_exon_len + $exon->end - $feature->end
	              + VC_MINUS_STRAND_HACK_BP;
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
  my $hit = $sorted_nuc_evidence_arr[0];
  my $prev_exon = $$hit{exon};
  $prev_hseqname = '-' x 1000;	# fake initial ID
  for (my $i = 0; $i < @sorted_nuc_evidence_arr; $i++) {
    my $hit = $sorted_nuc_evidence_arr[$i];
    my $hseq_str = $$hit{hseq};
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

  # change case at exon boundaries

  my @exon_lengths = ();
  foreach my $exon (@all_exons) {
    push @exon_lengths, ($exon->end - $exon->start + 1);
  }
  foreach my $evidence_line (@filtered_evidence_arr) {
    my $seq_str = $evidence_line->seq;
    $total_exon_len = 0;
    for (my $i = 0; $i < @exon_lengths; $i++) {
      my $replacement = substr $seq_str, $total_exon_len,
                               $exon_lengths[$i];
      if (! ($i & 2)) {
        $replacement = uc $replacement;
      } else {
        $replacement = lc $replacement;
      }
      substr $seq_str, $total_exon_len, $exon_lengths[$i], $replacement;
      $total_exon_len += $exon_lengths[$i];
    }
    $evidence_line->seq($seq_str);
  }

  return \@filtered_evidence_arr;
}

1;

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
                          -TRANSCRIPTID => $tr_stable_id,
			  -PFETCH       => '/path/to/pfetch');
 my $alignment_arr_ref = $ea->fetch_alignment;
 $ea->transcriptid($other_tr_dbID);
 my $other_alignment_arr_ref = $ea->fetch_alignment;
 my $short_alignment_arr_ref = $ea->fetch_alignment($hid);
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

# modify VCs' features' genomic start/end by the following:
use constant VC_HACK_BP  => +1;

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
                                -TRANSCRIPTID => $transcript_id,
				-PFETCH       => '/path/to/pfetch');
		my $cont_ea = Bio::EnsEMBL::EvidenceAlignment->new(
                                -DBADAPTOR => $dba,
                                -CONTIGID  => $contig_stable_id);
    Function:   Initialises EvidenceAlignment object
    Returns :   An EvidenceAlignment object
    Args    :   Database adaptor object and an ID string (-CONTIGID
                with contig ID or -TRANSCRIPTID with transcript
		ID, which may be either the transcript stable ID or
		the transcript dbID); optional full path of pfetch
		binary (-PFETCH), defaulting to whatever is in the
		user's search path.

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($transcriptid, $contigid, $dbadaptor, $pfetch) = $self->_rearrange(
                                                         ['TRANSCRIPTID',
						          'CONTIGID',
                                                          'DBADAPTOR',
						          'PFETCH'],
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
  if (! $pfetch) {
    $pfetch = 'pfetch';
  }
  $self->pfetch($pfetch);
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

=head2 pfetch

    Title   :   pfetch
    Usage   :   $ea->pfetch('/path/to/pfetch');
    Function:   get/set for pfetch, the location of the pfetch
                binary.

=cut

sub pfetch {
  my $obj = shift;
  if( @_ ) {
    my $value = shift;
    $obj->{evidencealignment_pfetch} = $value;
  }
  return $obj->{evidencealignment_pfetch};
}

=head2 transcriptid

    Title   :   transcriptid
    Usage   :   $ea->transcriptid($transcript_id);
    Function:   get/set for transcript id string, which may be
                either the stable ID or the internal (db) ID
		(setting transcriptid also unsets contigid)

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
    Function:   get/set for contig stable id string (setting
                contigid also unsets transcriptid)

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
                $ea->_get_features_from_transcript($transcript_obj, $vc, $hid);
    Function:   use SGP adaptor supplied to get evidence off a VC
                of the transcript supplied; features not overlapping
		any exon are cut; if a list of hit accession numbers
		are given, features not involving those accession
		numbers are cut; genomic start and end are always
		modified by VC_HACK_BP; duplicate features are
		removed
    Returns :   array of featurepairs

=cut

sub _get_features_from_transcript {
  my ($self, $transcript_obj, $vc) = splice @_, 0, 3;
  $self->throw('interface fault') if (!$self or !$transcript_obj or !$vc);
  my @wanted_arr = @_;

  my @exons = $transcript_obj->get_all_Exons;
  my $strand = $exons[0]->strand;
  my @all_features = $vc->get_all_SimilarityFeatures;
  
  my @features = ();
  FEATURE_LOOP:
  foreach my $feature (@all_features) {
    next unless $feature->primary_tag =~ /similarity/i;
    if ($feature->strand == $strand) {
      # fix VC-related coordinate problem
      $feature->start($feature->start + VC_HACK_BP);
      $feature->end($feature->end + VC_HACK_BP);
      # store on overlap, unless we're not interested in this hit seq.
      foreach my $exon (@exons) {
        if ($exon->overlaps($feature)) {
	  my $hid = $feature->hseqname;
	  my $wanted = 0;
	  if (@wanted_arr) {
	    foreach (@wanted_arr) {	# we only want some hits
	      if ($_ eq $hid) {
	        $wanted = 1;
	      }
	    }
	  } else {			# we want all hits
	    $wanted = 1;
	  }
	  if ($wanted) {
	    push @features, $feature;
	    next FEATURE_LOOP;
	  }
	}
      }
    }
  }

  # remove duplicates
  my @sorted_features = sort {    $a->hseqname cmp $b->hseqname
                               || $a->start    <=> $b->start
                               || $a->end      <=> $b->end
                               || $a->hstart   <=> $b->hstart
                               || $a->hend     <=> $b->hend
                               || $a->strand   <=> $b->strand
			     } @features;
  for (my $i = 1; $i < @sorted_features; $i++) {
    my $f1 = $sorted_features[$i];
    my $f2 = $sorted_features[$i-1];
    if ( not $f1->hseqname cmp $f2->hseqname
          || $f1->start    <=> $f2->start
          || $f1->end      <=> $f2->end
          || $f1->hstart   <=> $f2->hstart
          || $f1->hend     <=> $f2->hend
          || $f1->strand   <=> $f2->strand )
    {
      splice @sorted_features, $i, 1;
      $i--;
    }
  }
  
  return @sorted_features;
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
    next unless $feature->primary_tag =~ /similarity/i;
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

=head2 _get_Seqs_by_accs

    Title   :   _get_Seqs_by_accs
    Usage   :   $seqs = $self->_get_Seqs_by_accs(@hid_arr);
    Function:   Does the sequence retrieval for an array of
                accesions in one pfetch call. This is closely
		based on
		Bio/EnsEMBL::Pipeline::SeqFetcher::Pfetch in
		ensembl-pipeline. Having it here allows
		the current module to work for both pipeline
		and Web users.
    Returns :   Array of Bio::Seq

=cut

sub _get_Seqs_by_accs {
  my ($self, @acc) = @_;

  if (!defined(@acc) || scalar(@acc < 1)) {
    $self->throw("No accession input");
  }

  my @seq;
  my $newseq;
  my $tracker = 0;
  my $pfetch = $self->pfetch;

  my $command = "$pfetch -q ";
  $command .= join " ", @acc;

  open(EVIDENCEALIGNMENT_IN_FH,"$command |")
    or $self->throw("Error opening pipe to pfetch for $pfetch");
  while(<EVIDENCEALIGNMENT_IN_FH>){

    chomp;
    eval{
      if(defined $_ && $_ ne "no match") {
        $newseq = new Bio::Seq('-seq'               => $_,
                               '-accession_number'  => $acc[$tracker],
                               '-display_id'        => $acc[$tracker]);
      }
    };

    if($@){
      print STDERR "$@\n";
    }

    if (defined $newseq){
      push (@seq, $newseq);
    }
    else{
      $self->warn("Could not even pfetch sequence for [" . $acc[$tracker] . "]\n");
    }
    $tracker++;
  }

  close EVIDENCEALIGNMENT_IN_FH
    or $self->throw("Error running pfetch for $pfetch");
  return @seq;
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
                my $seq_arr_ref = $ea->fetch_alignment($hid);
    Function:   gets transcript or raw contig and corresponding
                evidence or similarity features; for raw
		contigs, these are displayed for the forward
		strand followed by the reverse strand
    Args    :   for transcripts, an optional list of accession
		numbers of hit sequences of interest; if none
		are given, all relevant hit sequences are
		retrieved; note that for contigs, all relevant
		hit sequences are always retrieved
    Returns :   reference to array of Bio::PrimarySeq

=cut

sub fetch_alignment {
  my ($self) = shift;
  $self->throw('interface fault') if (! $self);
  $self->throw('must have a transcript or contig ID and a DB adaptor object')
    unless (($self->transcriptid || $self->contigid) && $self->dbadaptor);

  if ($self->transcriptid) {
    return $self->_get_aligned_evidence_for_transcript($self->transcriptid,
                                                       $self->dbadaptor, @_);
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

  my $clump_size = 1000;	# number of sequences to fetch at once
  my %hits_hash = ();
  my @hseqnames = ();
  for (my $i = 0; $i < @$features_arr_ref; $i++) {
    my $hseqname = $$features_arr_ref[$i]->hseqname;
    if (! exists $hits_hash{$hseqname})
    {
      push @hseqnames, $hseqname;
      if ((@hseqnames % $clump_size) == 0) {
        my @tofetch = sort { $a cmp $b } @hseqnames;
	for (my $i = 1; $i < @tofetch; $i++) {
          if ($tofetch[$i] eq $tofetch[$i-1]) {
	    splice @tofetch, $i, 1;
	    $i--;
	  }
	}
	for (my $i = 1; $i < @tofetch; $i++) {
	  if ($tofetch[$i] =~ /\|/) {	# pipe would make trouble
	    splice @tofetch, $i, 1;
	    $i--;
	  }
	}
        my @seqs;
	if (@tofetch) {
	  @seqs = $self->_get_Seqs_by_accs(@tofetch);
	}
	foreach my $seq_obj (@seqs) {
          $hits_hash{$seq_obj->accession_number} = $seq_obj;
	}
	foreach my $hseqname (@hseqnames) {
	  if (not exists $hits_hash{$hseqname}) {
	    $self->warn("couldn't fetch sequence for hit $hseqname");
	    $hits_hash{$hseqname} = undef;
	  }
	}
        @hseqnames = ();
      }
    }
  }
  if (@hseqnames) {	# fetch the non-clump-sized remainder
    my @tofetch = sort { $a cmp $b } @hseqnames;
    for (my $i = 1; $i < @tofetch; $i++) {
      if ($tofetch[$i] eq $tofetch[$i-1]) {
        splice @tofetch, $i, 1;
        $i--;
      }
    }
    for (my $i = 1; $i < @tofetch; $i++) {
      if ($tofetch[$i] =~ /\|/) {	# pipe would make trouble
        splice @tofetch, $i, 1;
        $i--;
      }
    }
    my @seqs;
    if (@tofetch) {
      @seqs = $self->_get_Seqs_by_accs(@tofetch);
    }
    foreach my $seq_obj (@seqs) {
      $hits_hash{$seq_obj->accession_number} = $seq_obj;
    }
    foreach my $hseqname (@hseqnames) {
      if (not exists $hits_hash{$hseqname}) {
        $self->warn("couldn't fetch sequence for hit $hseqname");
        $hits_hash{$hseqname} = undef;
      }
    }
  }

  return \%hits_hash;
}

# _evidence_lines_sort: takes reference to an array of evidence lines
# and reference to a hash of scores,
# returns reference to the former sorted by score (highest score first),
# with ties sorted alphabetically

sub _evidence_lines_sort {
  my ($self, $tmp_evidence_arr_ref, $scores_hash_ref) = @_;
  $self->throw('interface fault') if (@_ != 3);

  my @sorted_arr = sort {
    $$scores_hash_ref{$b->accession_number}
      <=> $$scores_hash_ref{$a->accession_number}
      ||  $a->accession_number cmp $b->accession_number
  } @$tmp_evidence_arr_ref;

  return \@sorted_arr;
}

# _get_per_hid_effective_scores: takes reference to an array of features,
# returns a reference to a hash giving the 'effective score' for each
# hseqname, i.e., the score-like value upon which we wish to sort.

sub _get_per_hid_effective_scores {
  my ($self, $feature_arr_ref) = @_;
  $self->throw('interface fault') if (@_ != 2);

  my %per_hid_effective_scores = ();
  foreach my $feature (@$feature_arr_ref) {
    my $hseqname = $feature->hseqname;
    my $feature_len = $feature->end - $feature->start + 1;
    if (not exists $per_hid_effective_scores{$hseqname})
    {
      $per_hid_effective_scores{$hseqname} = 0;
    }
    $per_hid_effective_scores{$hseqname} += $feature_len;
    # another possibility:
    # $per_hid_effective_scores{$hseqname} += $feature_len * $feature->score
  }
  return \%per_hid_effective_scores;
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
  my $per_hid_effective_scores_hash_ref
    = $self->_get_per_hid_effective_scores(\@features);
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
    if (($feature->start < 1) || ($feature->end > $dna_len_bp)) {
      $self->warn("genomic coordinates out of range: start " .
        $feature->start . ", end " . $feature->end);
      next PEP_FEATURE_LOOP;
    }
    my $hit_seq_obj = $$hits_hash_ref{$feature->hseqname};
    if (! $hit_seq_obj) {
      next PEP_FEATURE_LOOP;	# already warned in _get_hits()
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
  my @tmp_pep_evidence_arr = ();

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
      push @tmp_pep_evidence_arr, $evidence_obj;
    }
    $prev_hseqname = $$hit{hseqname};
  }

  my $sorted_pep_evidence_lines_ref =
    $self->_evidence_lines_sort(\@tmp_pep_evidence_arr,
                                $per_hid_effective_scores_hash_ref);
  push @evidence_arr, @$sorted_pep_evidence_lines_ref;

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
    if (($feature->start < 1) || ($feature->end > $dna_len_bp)) {
      $self->warn("genomic coordinates out of range: start " .
        $feature->start . ", end " . $feature->end);
      next NUC_FEATURE_LOOP;
    }
    my $hit_seq_obj = $$hits_hash_ref{$feature->hseqname};
    if (! $hit_seq_obj) {
      next NUC_FEATURE_LOOP;	# already warned in _get_hits()
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
  my @tmp_nuc_evidence_arr = ();

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
      push @tmp_nuc_evidence_arr, $evidence_obj;
    }
    $prev_hseqname = $$hit{hseqname};
  }

  my $sorted_nuc_evidence_lines_ref =
    $self->_evidence_lines_sort(\@tmp_nuc_evidence_arr,
                                $per_hid_effective_scores_hash_ref);
  push @evidence_arr, @$sorted_nuc_evidence_lines_ref;

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
# a DB adaptor, and an optional list of hit sequence accession numbers
# returns ref to an array of Bio::PrimarySeq

sub _get_aligned_evidence_for_transcript {
  my ($self, $transcript_id, $db) = splice @_, 0, 3;
  $self->throw('interface fault') if (!$self or !$transcript_id or !$db);

  my $sgp = $db->get_StaticGoldenPathAdaptor;
  my @evidence_arr;	# a reference to this is returned
  my $evidence_obj;
  my $ta = $db->get_TranscriptAdaptor;
  my $ga = $db->get_GeneAdaptor;
  my $ea = $db->get_ExonAdaptor;
  my $transcript_name_to_display = $transcript_id;

  # get all exons in VC coordinates
  my $transcript_obj;
  my $transcript_dbID;
  if ($transcript_id =~ /^ENS/i) {	# stable ID
    $transcript_obj = $ta->fetch_by_stable_id($transcript_id);
    $transcript_dbID = $transcript_obj->dbID;
  } else {	# internal (db) ID
    $transcript_dbID = $transcript_id;
    $transcript_obj = $ta->fetch_by_dbID($transcript_dbID);
  }
  my $vc = $sgp->fetch_VirtualContig_of_transcript_by_dbID($transcript_dbID,
                                                          1000);
  my @all_exons = $transcript_obj->get_all_Exons;
  foreach my $exon (@all_exons) {
    $vc->_convert_seqfeature_to_vc_coords($exon);
  }

  my @features = $self->_get_features_from_transcript($transcript_obj, $vc, @_);
  my $per_hid_effective_scores_hash_ref =
    $self->_get_per_hid_effective_scores(\@features);
  my $hits_hash_ref = $self->_get_hits(\@features);
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

  # get 3' UTR length
  my $total_3prime_utr_len =   $cdna_len_bp
                             - length($transcript_obj->translate->seq) * 3
			     - $total_5prime_utr_len;
  if ($total_3prime_utr_len < 0) {
    $total_3prime_utr_len = 0;	# disaster recovery
  }

  # get translation, avoiding transcript's translate method due to rare problem
  my $seq_to_translate = substr $nucseq_str, $total_5prime_utr_len;
  my $cdna_obj = Bio::PrimarySeq->new(
                    -seq              => $seq_to_translate,
                    -id               => 0,
                    -accession_number => $transcript_name_to_display,
		    -moltype          => 'dna'
		  );
  my $translation_including_3prime_utr = $cdna_obj->translate->seq;
  
  # translation itself forms our first row of 'evidence'
  my $evidence_line = $translation_including_3prime_utr;
  $evidence_line = $self->_pad_pep_str($evidence_line);
  # adjust evidence line for UTRs
  $evidence_line = ('-' x $total_5prime_utr_len) . $evidence_line;
  $evidence_line = substr $evidence_line, 0,   $cdna_len_bp
                                             - $total_3prime_utr_len;
  $evidence_line .= '-' x $total_3prime_utr_len;
  # store
  $evidence_obj = Bio::PrimarySeq->new(
                    -seq              => $evidence_line,
                    -id               => 0,
                    -accession_number => $transcript_name_to_display,
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
        unless $exon->overlaps($feature);
      my $hit_seq_obj = $$hits_hash_ref{$feature->hseqname};
      if (! $hit_seq_obj) {
        next PEP_FEATURE_LOOP;	# already warned in _get_hits()
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
      if ($feature->start < $exon->start) {
        my $old_flen = $flen;
	$flen -= $exon->start - $feature->start;
        $feature->start($exon->start);
	if ($exon->strand > 0) {	# trim start of hit
          $hseq = substr $hseq, $flen - $old_flen, $flen;
	} else {			# trim end of hit
          $hseq = substr $hseq, 0, $flen;
	}
      }
      if ($feature->end > $exon->end) {
        my $old_flen = $flen;
	$flen -= $feature->end - $exon->end;
	$feature->end($exon->end);
        if ($exon->strand > 0) {	# trim end of hit
          $hseq = substr $hseq, 0, $flen;
	} else {			# trim start of hit
	  $hseq = substr $hseq, $old_flen - $flen, $flen;
	}
      }
      my $hindent_bp;
      if ($exon->strand > 0) {
        $hindent_bp = $total_exon_len + $feature->start - $exon->start;
      } else {
        $hindent_bp = $total_exon_len + $exon->end - $feature->end;
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
  my @tmp_pep_evidence_arr = ();

  $evidence_line = '';
  my $prev_hseqname = '-' x 1000;	# fake initial ID
  for (my $i = 0; $i < @sorted_pep_evidence_arr; $i++) {
    my $hit = $sorted_pep_evidence_arr[$i];
    if ($$hit{hseqname} ne $prev_hseqname) {	# make new evidence line
      $evidence_line = '-' x $cdna_len_bp;
    }

    # splice in the evidence fragment
    my $hseqlen = length $$hit{hseq};
    substr $evidence_line, $$hit{hindent}, $hseqlen, $$hit{hseq};

    # store if end of evidence line
    if (($i == $#sorted_pep_evidence_arr)
     || ($sorted_pep_evidence_arr[$i+1]{hseqname} ne $$hit{hseqname}))
    {
      # purge the UTRs of protein 'evidence'
      substr $evidence_line, 0, $total_5prime_utr_len,
             ('-' x $total_5prime_utr_len);
      substr $evidence_line, $cdna_len_bp - $total_3prime_utr_len,
             $total_3prime_utr_len, ('-' x $total_3prime_utr_len);
      $evidence_obj = Bio::PrimarySeq->new(
                      -seq              => $evidence_line,
                      -id               => 0,
  		      -accession_number => $$hit{hseqname},
		      -moltype          => $$hit{moltype}
		    );
      push @tmp_pep_evidence_arr, $evidence_obj;
    }
    $prev_hseqname = $$hit{hseqname};
  }

  my $sorted_pep_evidence_lines_ref =
    $self->_evidence_lines_sort(\@tmp_pep_evidence_arr,
                                $per_hid_effective_scores_hash_ref);
  push @evidence_arr, @$sorted_pep_evidence_lines_ref;

  # nucleic acid evidence

  # cDNA itself forms first row of nucleic acid evidence
  $evidence_obj = Bio::PrimarySeq->new(
                    -seq              => $nucseq_str,
                    -id               => 0,
     		    -accession_number => $transcript_name_to_display,
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
        unless $exon->overlaps($feature);
      my $hit_seq_obj = $$hits_hash_ref{$feature->hseqname};
      if (! $hit_seq_obj) {
	next NUC_FEATURE_LOOP;	# already warned in _get_hits()
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
      if ($feature->start < $exon->start) {
        my $old_flen = $flen;
	$flen -= $exon->start - $feature->start;
	$hlen -= $exon->start - $feature->start;
        $feature->start($exon->start);
	if ($exon->strand > 0) {	# trim start of hit
	  $feature->hstart($feature->hstart + $old_flen - $flen);
	} else {			# trim end of hit
	  $feature->hend($feature->hend - ($old_flen - $flen));
	}
      }
      if ($feature->end > $exon->end) {
        my $old_flen = $flen;
	$flen -= $feature->end - $exon->end;
	$hlen -= $feature->end - $exon->end;
	$feature->end($exon->end);
        if ($exon->strand > 0) {	# trim end of hit
	  $feature->hend($feature->hend - ($old_flen - $flen));
	} else {			# trim start of hit
	  $feature->hstart($feature->hstart + $old_flen - $flen);
	}
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
        $hindent_bp =   $total_exon_len + $feature->start - $exon->start;
      } else{
        $hindent_bp =    $total_exon_len + $exon->end - $feature->end;
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
  my @tmp_nuc_evidence_arr = ();

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
      push @tmp_nuc_evidence_arr, $evidence_obj;
    }
    $prev_hseqname = $$hit{hseqname};
    $prev_exon = $$hit{exon};
  }

  my $sorted_nuc_evidence_lines_ref =
    $self->_evidence_lines_sort(\@tmp_nuc_evidence_arr,
                                $per_hid_effective_scores_hash_ref);
  push @evidence_arr, @$sorted_nuc_evidence_lines_ref;

  # remove blank evidence lines

  my @filtered_evidence_arr = ();
  foreach my $evidence_line (@evidence_arr) {
    push @filtered_evidence_arr, $evidence_line
      if ($$evidence_line{seq} =~ /[^-]/);
  }

  # for transcripts: remove cDNA if one protein hid specified,
  # remove translation if one nucleotide hid specified, but (in
  # desperation) leave translation in place if there are no hits
  # of any kind

  if ($self->transcriptid) {
    my $prot_lines = 0;
    my $nuc_lines = 0;
    foreach my $evidence_line (@filtered_evidence_arr) {
      if ($evidence_line->moltype eq 'protein') {
        $prot_lines++;
      } else {
        $nuc_lines++;
      }
    }
    if ($nuc_lines eq 1) {	# for nucleotides, we have cDNA only
      splice @filtered_evidence_arr, $#filtered_evidence_arr, 1;
    }
    if (@filtered_evidence_arr > 1 and $prot_lines == 1) {
      splice @filtered_evidence_arr, 0, 1;
    }
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
      if ($i % 2) {
        $replacement = lc $replacement;
      } else {
        $replacement = uc $replacement;
      }
      substr $seq_str, $total_exon_len, $exon_lengths[$i], $replacement;
      $total_exon_len += $exon_lengths[$i];
    }
    $evidence_line->seq($seq_str);
  }

  return \@filtered_evidence_arr;
}

1;

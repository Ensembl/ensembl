# EnsEMBL RawContig object
#
# Copyright EMBL-EBI 2001
#
# cared for by:: Arne Stabenau
# Date : 04.12.2001
#


=head1 NAME

Bio::EnsEMBL::RawContig
  Contig object which represents part of an EMBL Clone.Mostly for 
  database usage

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut


package Bio::EnsEMBL::RawContig;

use vars qw( @ISA );
use strict;


use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::PrimarySeq;

@ISA = qw( Bio::EnsEMBL::Root Bio::PrimarySeqI );

=head2 new

  Arg  1    : int dbID
  Arg  2    : Bio::EnsEMBL::DBSQL::RawContigAdaptor adaptor
  Arg  3    : txt contigName
  Arg  4    : Bio::SeqI sequenceContainer
  Arg  5    : int sequenceLength
  Arg  6    : Bio::EnsEMBL::Clone clone
  Arg  7    : int emblCloneBaseOffset
  Arg  7    : int emblCloneContigNumber
  Arg  8    : int emblCloneBaseOffset
  Function  : creates RawContig. Neds either dbID and Adaptor or clone and sequence.
              With dbID its connected to DB, with the other it may be stored. 
  Returntype: Bio::EnsEMBL::RawContig
  Exceptions: none
  Caller    : RawContigAdaptor, Data upload script

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = {};
  bless $self, $class;
  
  my ( $dbID, $adaptor, $name, $sequence, $length,
       $clone, $offset ) = @args;


  (defined $dbID) && $self->dbID( $dbID );
  (defined $adaptor) && $self->adaptor( $adaptor );
  (defined $clone) && $self->clone( $clone );
  (defined $sequence) && $self->sequence( $sequence );
  (defined $name) && $self->name( $name );
  (defined $length) && $self->length( $length );
  (defined $offset) && $self->embl_offset( $offset );

  return $self;
}



sub adaptor {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_adaptor} = $arg );

  return $self->{_adaptor};
}



sub dbID {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_dbID} = $arg );
  
  return $self->{_dbID};
}

    

sub name {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_name} = $arg ;
  } else {
    if( ! defined $self->{_name} &&
      defined $self->adaptor() ) {
      $self->adaptor->fetch_attributes( $self );
    }
  }
  
  return $self->{_name};
}


# scp: RawContig->primary_seq is used by the pipeline/bioperl
# (which expects an object to be returned that implements 'id'.
# As RawContig->primary_seq now returns a RawContig, this
# object needs to implement 'id' also.

sub id {
  my ($self) = shift;

  return $self->name || $self->dbID;
}



sub length {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_length} = $arg ;
  } else {
    if( ! defined $self->{_length} &&
	defined $self->adaptor() ) {
      $self->adaptor->fetch_attributes( $self );
    }
  }
  
  return $self->{_length};
}



=head2 seq

  Arg [1]    : none
  Example    : $dna = $contig->seq();
  Description: Obtains the dna sequence of this entire contig as a string
               (positive strand).
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub seq {
  my ($self, $arg) = @_;

  # Sequence can be set manually
  if($arg) {
    $self->{_seq} = $arg;
  }
  if($self->{_seq}) {
    return $self->{_seq};
  }

  #or retrieved from the database
  if($self->adaptor()) {
    my $sa = $self->adaptor->db->dnadb->get_SequenceAdaptor(); 
    return $sa->fetch_by_RawContig_start_end_strand($self, 1, -1, 1);
  }
  
  $self->warn("RawContig seq not set, and no db is available");
  return '';
}



=head2 subseq

  Arg [1]    : int $start
               The start basepair of the sequence to obtain
  Arg [2]    : int $end
               The end basepair of the sequence to obtain
  Arg [3]    : (optional) int $strand
               The strand of the sequence to obtain.  Default is the contigs
               forward strand.
  Example    : $dna = $contig->subseq(1, 1000, -1);
  Description: Obtains a subsequence of this contigs dna sequence
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub subseq {
  my ($self, $start, $end, $strand) = @_;

  if ( $end < $start ) {
    $self->throw("End coord is less then start coord to call on RawContig " .
		 "subseq.");
  }

  if ( !defined $strand || ( $strand != -1 && $strand != 1 )) {
    $strand = 1;
  }

  #if the sequence of this contig has been manually set retrieve its substring
  if(my $str = $self->{_seq}) {
    $str = substr($str, $start -1, $end - $start + 1);

    if($strand == -1) {
      $str = reverse $str;
      $str =~ tr/ACTGactg/TGACtgac/;
    }

    return $str;
  }
   
  unless($self->adaptor) {
    $self->warn("RawContig::subseq no sequence set and no db available");
    return '';
  }

  my $seq_adaptor = $self->adaptor->db->get_SequenceAdaptor();

  my $sub_seq = $seq_adaptor->fetch_by_RawContig_start_end_strand($self, 
							       $start, $end, 
							       $strand);
  return $sub_seq;
}



=head2 get_base_count

  Arg [1]    : none
  Example    : $gc_content = $contig->get_base_count()->{'%gc'};
  Description: Retrieves a hashref containing the counts of each bases in the
               sequence spanned by this slice.  The format of the hash is :
               { 'a' => num,
                 'c' => num,
                 't' => num,
                 'g' => num,
                 'n' => num,
                 '%gc' => num }
               
               All bases which are not in the set [A,a,C,c,T,t,G,g] are 
               included in the 'n' count.  The 'n' count could therefore be
               inclusive of ambiguity codes such as 'y'.
               The %gc is the ratio of GC to AT content as in:
               total(GC)/total(ACTG) * 100
  Returntype : hashref
  Exceptions : none
  Caller     : general

=cut

sub get_base_count {
  my $self = shift;

  my $seq = $self->seq();
  my $len = $self->length();

  my $a = $seq =~ tr/Aa/Aa/;
  my $c = $seq =~ tr/Cc/Cc/;
  my $t = $seq =~ tr/Tt/Tt/;
  my $g = $seq =~ tr/Gg/Gg/;

  my $gc_content = 0;
  if($a || $g || $c || $t) {  #avoid divide by 0
    $gc_content = sprintf( "%1.2f", (($g + $c)/($a + $g + $t + $c)) * 100);
  }

  return {'a' => $a,
	  'c' => $c,
	  't' => $t,
	  'g' => $g,
	  'n' => $len - $a - $c - $t - $g,
	  '%gc' => $gc_content};
}


=head2 get_repeatmasked_seq

  Arg [1]    : string \@logic_names (optional)
  Arg [2]    : int $soft_masking_enable (optional)
  Example    : $slice->get_repeatmasked_seq or
               $slice->get_repeatmasked_seq(['RepeatMask'],1)
  Description: Returns Bio::PrimarySeq containing the masked 
               (repeat replaced by N) 
               or soft-masked (when Arg[2]=1, repeat in lower case while 
               non repeat in upper case) sequence corresponding to the 
               Slice object. Will only work with database connection to get 
               repeat features.
  Returntype : Bio::PrimarySeq
  Exceptions : none
  Caller     : general.

=cut

sub get_repeatmasked_seq {
    my ($self, $logic_names, $soft_mask) = @_;
    
    unless ($logic_names && @$logic_names) {
        $logic_names = [ '' ];
    }

    unless (defined $soft_mask) {
      $soft_mask = 0;
    }

    my $repeats = [];

    foreach my $l (@$logic_names) {
	push @{$repeats}, @{$self->get_all_RepeatFeatures($l)};
    }

    my $dna = $self->seq();
    my $masked_dna = $self->_mask_features($dna,$repeats,$soft_mask);
    my $masked_seq = Bio::PrimarySeq->new('-seq'        => $masked_dna,
					  '-display_id' => $self->name,
					  '-primary_id' => $self->name,
					  '-moltype'    => 'dna'
					 );
    return $masked_seq;
}


=head2 _mask_features

  Arg [1]    : string $dna_string
  Arg [2]    : array_ref \@repeats
               reference to a list Bio::EnsEMBL::RepeatFeature
               give the list of coordinates to replace with N or with lower case
  Arg [3]    : int $soft_masking_enable (optional)
  Example    : 
  Description: replaces string positions described in the RepeatFeatures
               with Ns (default setting), or with the lower case equivalent (soft masking)
  Returntype : string 
  Exceptions : none
  Caller     : get_repeatmasked_seq

=cut

sub _mask_features {
  my ($self,$dnastr,$repeats,$soft_mask) = @_;
    
  # explicit CORE::length call, to avoid confusion with Slice::length method
  my $dnalen = CORE::length($dnastr);
    
 REP:foreach my $f (@{$repeats}) {
      
    my $start  = $f->start;
    my $end    = $f->end;
    my $length = ($end - $start) + 1;
    
    # check if we get repeat completely outside of expected slice range
    if ($end < 1 || $start > $dnalen) {
      $self->warn("Repeat completely outside RawContig coordinates!".
		  "That should not happen! repeat_start $start or " .
		  "repeat_end $end not within [1-$dnalen] RawContig " .
		  "range coordinates\n");
      next REP;
    }
    
    # repeat partly outside slice range, so correct
    # the repeat start and length to the slice size if needed
    if ($start < 1) { 
      $start = 1;
      $length = ($end - $start) + 1;
    }

    # repeat partly outside slice range, so correct
    # the repeat end and length to the slice size if needed
    if ($end > $dnalen) {
      $end = $dnalen;
      $length = ($end - $start) + 1;
    }
    
    $start--;
    
    my $padstr;
    
    if ($soft_mask) {
      $padstr = lc substr ($dnastr,$start,$length);
    } else {
      $padstr = 'N' x $length;
    }
    substr ($dnastr,$start,$length) = $padstr;
  }
  return $dnastr;
}


=head2 get_all_PredictionTranscripts

  Args      : none
  Function  : connect to database through set adaptor and retrieve the 
              PredictionFeatures for this contig.
  Returntype: listref Bio::EnsEMBL::PredictionTranscript 
              (previously this returned a SeqFeature)
  Exceptions: none
  Caller    : general

=cut

sub get_all_PredictionTranscripts {
  my $self = shift;
  my $logic_name = shift;
  
  if( ! defined $self->adaptor() ) {
    $self->warn( "Need db connection for get_all_PredictionFeatures()" );
    return ();
  }
  
  my $pta = $self->adaptor->db->get_PredictionTranscriptAdaptor();
    
  return $pta->fetch_all_by_RawContig($self, $logic_name);
}



=head2 get_all_RepeatFeatures

  Args      : none
  Function  : connect to database through set adaptor and retrieve the 
              repeatfeatures for this contig.
  Returntype: listref Bio::EnsEMBL::RepeatFeature
  Exceptions: none
  Caller    : general, get_repeatmasked_seq()

=cut

sub get_all_RepeatFeatures {
   my $self = shift;
   my $logic_name = shift;
  
   if( ! defined $self->adaptor() ) {
     $self->warn( "Need db connection for get_all_RepeatFeatures()" );
     return ();
   }

   my $rfa = $self->adaptor()->db->get_RepeatFeatureAdaptor();

   return $rfa->fetch_all_by_RawContig( $self , $logic_name);
}



=head2 get_all_SimilarityFeatures

  Args      : none
  Function  : connect to database through set adaptor and retrieve the 
              SimilarityFeatures for this contig.
  Returntype: list Bio::EnsEMBL::FeaturePair
  Exceptions: none
  Caller    : general

=cut

sub get_all_SimilarityFeatures {
  my ($self, $logic_name, $score) = @_;
  
  if( ! defined $self->adaptor() ) {
    $self->warn( "Need db connection for get_all_SimilarityFeatures()" );
    return ();
  }
  
  my @out;
  my $dafa = $self->adaptor->db->get_DnaAlignFeatureAdaptor();
  my $pafa = $self->adaptor->db->get_ProteinAlignFeatureAdaptor();
  push @out, @{$dafa->fetch_all_by_RawContig_and_score($self, $score, $logic_name)};
  push @out, @{$pafa->fetch_all_by_RawContig_and_score($self, $score, $logic_name)};
    
  return \@out;
}


=head2 get_all_DnaAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the dna align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Example    : @dna_align_feats = $contig->get_all_DnaAlignFeatures()
  Description: Retrieves the DnaDnaAlignFeatures which overlap this contig
  Returntype : list of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : none
  Caller     :general

=cut

sub get_all_DnaAlignFeatures {
   my ($self, $logic_name, $score) = @_;


   if( ! defined $self->adaptor() ) {
     $self->warn( "Need db connection for get_all_DnaAlignFeatures()" );
     return ();
   }

   my $dafa = $self->adaptor->db->get_DnaAlignFeatureAdaptor();

   return $dafa->fetch_all_by_RawContig_and_score($self,$score, $logic_name);
}



=head2 get_all_ProteinAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the protein align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Example    : @pep_align_feats = $contig->get_all_ProteinAlignFeatures()
  Description: Retrieves the PepDnaAlignFeatures which overlap this contig.
  Returntype : list of Bio::EnsEMBL::PepDnaAlignFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_ProteinAlignFeatures {
  my ($self, $logic_name, $score) = @_;


  if( ! defined $self->adaptor() ) {
    $self->warn( "Need db connection for get_all_ProteinAlignFeatures()" );
    return ();
  }

  my $pafa = $self->adaptor()->db()->get_ProteinAlignFeatureAdaptor();

  return $pafa->fetch_all_by_RawContig_and_score($self, $score, $logic_name);
}



=head2 get_all_SimpleFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the simple features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Example    : @simple_feats = $contig->get_all_SimpleFeatures()
  Description: Retrieves the SimpleFeatures which overlap this contig.
  Returntype : list of Bio::EnsEMBL::DnaDnaAlignFeature
  Exceptions : none
  Caller     : general

=cut

sub get_all_SimpleFeatures {
  my ($self, $logic_name, $score) = @_;


  if( ! defined $self->adaptor() ) {
    $self->warn( "Need db connection for get_all_SimpleFeatures()" );
    return ();
  }

  my $sfa = $self->adaptor()->db()->get_SimpleFeatureAdaptor();

  return $sfa->fetch_all_by_RawContig_and_score($self, $score, $logic_name);
}


sub embl_offset {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_embl_offset} = $arg ;
  } else {
    if( ! defined $self->{_embl_offset} &&
      defined $self->adaptor() ) {
      $self->adaptor->fetch_attributes( $self );
    }
  }
  
  return $self->{_embl_offset};
}


sub clone {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_clone} = $arg ;
  } else {
    if( ! defined $self->{_clone} &&
      defined $self->adaptor() ) {
      if( !defined $self->_clone_id() ) {
	$self->adaptor->fetch_attributes($self);
      }
      $self->{_clone} = 
       $self->adaptor->db->get_CloneAdaptor->fetch_by_dbID($self->_clone_id());
    }
  }
  
  return $self->{_clone};
}

sub _clone_id {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_clone_id} = $arg ;
  }

  return $self->{_clone_id};
}




=head2 get_all_ExternalFeatures

  Arg [1]    : (optional) string $track_name
               If specified only features from ExternalFeatureAdaptors with 
               the track name $track_name are retrieved.  If not set, all 
               features from every ExternalFeatureAdaptor are retrieved.
  Example    : @xfeatures = @{$contig->get_all_ExternalFeatures};
  Description: Retrieves features on this contig from external feature adaptors
  Returntype : listref of Bio::SeqFeatureI implementing objects in contig
               coordinates 
  Exceptions : none
  Caller     : general

=cut

sub get_all_ExternalFeatures {
   my ($self, $track_name) = @_;

   my $features = [];

   my $xfa_hash = $self->adaptor->db->get_ExternalFeatureAdaptors;
   my @xf_adaptors = ();

   if($track_name) {
     #use a specific adaptor
     push @xf_adaptors, $xfa_hash->{$track_name};
   } else {
     #use all of the adaptors
     push @xf_adaptors, values %$xfa_hash;
   }

   foreach my $xfa (@xf_adaptors) {
     push @$features, @{$xfa->fetch_all_by_RawContig($self)};
   }

   return $features;
}



=head2 Methods included only for BioPerl compliance
=cut
###############################################################################

=head2 display_id

  Arg [1]    : none
  Example    : none
  Description: Only for BioPerl compliance.
  Returntype : string
  Exceptions : none
  Caller     : none

=cut

sub display_id{
  my $self = shift;

  return $self->id();
}

=head2 desc

  Arg [1]    : none
  Example    : none
  Description: Only for BioPerl compliance
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub desc{
  my $self = shift;
  return "RawContig, no description";
}

=head2 moltype

  Arg [1]    : none
  Example    : none
  Description: Only for BioPerl compliance
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub moltype {
  my $self = shift;
  return 'dna';
}

=head2 accession_number

  Arg [1]    : none
  Example    : none
  Description: Only for BioPerl compliance
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub accession_number {
  my $self = shift;
  return $self->dbID();
}



1;

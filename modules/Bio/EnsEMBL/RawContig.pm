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
use Bio::PrimarySeqI;

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

sub display_id {
  my $self = shift;
  $self->id();
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
    return $sa->fetch_by_Contig_start_end_strand($self, 1, -1, 1);
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
    $str = substr($str, $start -1, $end -1);

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

  my $seq_adaptor = $self->adaptor->db->dnadb->get_SequenceAdaptor();

  my $sub_seq = $seq_adaptor->fetch_by_Contig_start_end_strand($self, 
							       $start, $end, 
							       $strand);
  return $sub_seq;
}



=head2 get_repeatmasked_seq

  Args      : none
  Function  : Gives back the in memory repeatmasked seq. Will only work with
              Database connection to get repeat features.
  Returntype: Bio::PrimarySeq
  Exceptions: none
  Caller    : RunnableDB::Genscan::fetch_input(), other Pipeline modules.

=cut

sub get_repeatmasked_seq {
    my ($self, $logic_name) = @_;
    
    if(!$logic_name){
      $logic_name = 'RepeatMask';
    }

    my @repeats = $self->get_all_RepeatFeatures($logic_name);

    my $dna = $self->seq();
    my $masked_dna = $self->_mask_features($dna, @repeats);
    my $masked_seq = Bio::PrimarySeq->new(   '-seq'        => $masked_dna,
                                             '-display_id' => $self->name,
                                             '-primary_id' => $self->name,
                                             '-moltype' => 'dna',
					     );
    return $masked_seq;
}


=head2 _mask_features

  Arg  1    : txt $dna_string
  Arg  2    : list Bio::EnsEMBL::RepeatFeature @repeats
              list of coordinates to replace with N
  Function  : replaces string positions described im the RepeatFeatures
              with Ns. 
  Returntype: txt
  Exceptions: none
  Caller    : get_repeatmasked_seq

=cut

sub _mask_features {
    my ($self, $dnastr,@repeats) = @_;
    my $dnalen = CORE::length($dnastr);
    
  REP:foreach my $f (@repeats) {
      
      my $start    = $f->start;
      my $end	   = $f->end;
      my $length = ($end - $start) + 1;
      
      if ($start < 0 || $start > $dnalen || $end < 0 || $end > $dnalen) {
	  $self->warn("Coordinate mismatch - $start or $end not " .
	    "within $dnalen\n");
	  next REP;
      }
      
      $start--;
      
      my $padstr = 'N' x $length;
      
      substr ($dnastr,$start,$length) = $padstr;
  }
    return $dnastr;
}



=head2 get_all_RepeatFeatures

  Args      : none
  Function  : connect to database through set adaptor and retrieve the 
              repeatfeatures for this contig.
  Returntype: list Bio::EnsEMBL::RepeatFeature
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
   my @repeats = $rfa->fetch_by_Contig( $self , $logic_name);

   return @repeats;
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
  my ($self, $logic_name) = @_;
  
  if( ! defined $self->adaptor() ) {
    $self->warn( "Need db connection for get_all_SimilarityFeatures()" );
    return ();
  }
  
  my @out;
  my $dafa = $self->adaptor->db->get_DnaAlignFeatureAdaptor();
  my $pafa = $self->adaptor->db->get_ProteinAlignFeatureAdaptor();
      

  my @dnaalign = $dafa->fetch_by_Contig($self, $logic_name);
  my @pepalign = $pafa->fetch_by_Contig($self, $logic_name);
    
  push(@out, @dnaalign);
  push(@out, @pepalign);

  return @out;
}


=head2 get_all_PredictionFeatures

  Args      : none
  Function  : connect to database through set adaptor and retrieve the 
              PredictionFeatures for this contig.
  Returntype: list Bio::EnsEMBL::PredictionTranscript 
              (previously this returned a SeqFeature)
  Exceptions: none
  Caller    : general

=cut

sub get_all_PredictionFeatures {
  my $self = shift;
  
  if( ! defined $self->adaptor() ) {
    $self->warn( "Need db connection for get_all_PredictionFeatures()" );
    return ();
  }
  
  my $pta = $self->adaptor->db->get_PredictionTranscriptAdaptor();
    
  my @pred_feat = $pta->fetch_by_Contig($self);

  return @pred_feat;
}

sub accession_number {
  my $self = shift;
  
  $self->dbID();
}

sub moltype {
  return "DNA";
}

sub desc {
  return "Contig, no description";
}


=head2 get_all_ExternalFeatures

  Arg [1]    : none
  Example    : @external = $contig->get_all_ExternalFeatures
  Description: retrieves features generated by external feature factories
               attached to this database which are on this Contig.  
               See Bio::EnsEMBL::DB::ExternalFeatureFactoryI for details.
  Returntype : list of Bio::SeqFeatureI implementing objects 
  Exceptions : none
  Caller     : external

=cut

sub get_all_ExternalFeatures {
   my ($self) = @_;

   my @out = ();
   my $acc = $self->name();
   
#   $acc = $self->clone->id;

   my $offset = $self->embl_offset();

   unless($self->adaptor->db->can('_each_ExternalFeatureFactory')) {
     $self->warn("Database Adaptor is not capable of obtaining" .
		 "external features");
     return;
   }

   foreach my $extf ( $self->adaptor->db->_each_ExternalFeatureFactory ) {

     if( $extf->can('get_Ensembl_SeqFeatures_contig') ) {
       my @sfs = 
	 $extf->get_Ensembl_SeqFeatures_contig($self->dbID,
					       $self->clone->embl_version,
					       1,
					       $self->length,
		 			       $acc);
       push(@out,@sfs);
     }
     if( $extf->can('get_Ensembl_SeqFeatures_clone') ) {
       my @sfs = 
	 $extf->get_Ensembl_SeqFeatures_clone($acc, 
					      $self->clone->embl_version,
					      $self->embl_offset, 
					      $self->embl_offset+
					      $self->length());
       
       foreach my $sf (@sfs) {
	 my $start = $sf->start - $offset+1;
	 my $end   = $sf->end   - $offset+1;
	 if($start < 0 || $end > $self->length()) {
	   #discard features which are not entirely on this contig.
	   #features that span contigs should be shortened, but it is
           #unlikely anyone uses features spanning contigs anyway
	   next;
	 }
	 $sf->start($start);
	 $sf->end($end);
	 push(@out,$sf);
       }
     }
   }
   
   foreach my $f ( @out ) {
     $f->seqname($acc);
   }
   
   return @out;
}


=head2 sequence

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use seq method instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub sequence {
    my ($self,$arg) = @_;

    $self->warn("call to deprecated method Bio::EnsEMBL::RawContig::sequence "
		. "use the seq method instead\n");
 
    return $self->seq($arg);
}

=head2 dbobj

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use adaptor->db instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub dbobj {
   my ($self,$arg) = @_;

   $self->throw("RawContig::dbobj() is deprecated");

   return undef;
}


=head2 get_genscan_peptides

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_PredictionFeatures instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_genscan_peptides {
  my ($self, @args)  = @_;
  
  $self->warn("Use of deprecated method " .
	      "Bio::EnsEMBL::RawContig::get_genscan_peptides. " .
	      "Use get_PreedictionFeatures instead\n" );

  return $self->get_PredictionFeatures(@args);
}


=head2 primary_seq

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED Bio::EnsEMBL::RawContig is now a primary seq object.
              It should be used directly, or the seq method can be used
              to obtain a sequence string.
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub primary_seq {
  my $self = shift;

  $self->warn("Bio::EnsEMBL::RawContig is now a primary seq object use it " .
	      "directly or use the seq method to return a sequence string\n");

  return $self;
}


=head2 seq_old

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED should not have to use this method
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub seq_old {
  my $self = shift;
  my $arg = shift;

  $self->warn("call to deprecated method Bio::EnsEMBL::RawContig::seq_old");

  if( defined $arg ) {
    $self->{_seq} = $arg ;
  } else {
    if( ! defined $self->{_seq} &&
      defined $self->adaptor() ) {
	 print STDERR "Fetching sequence\n";
      $self->adaptor->fetch_attributes( $self );
    }
  }

  return $self->{_seq};
}


1;

# EnsEMBL Exon reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 22.11.2001
#

=head1 NAME

Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor - 
  MySQL Database queries to load and store PredictionExons

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut



package Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor;

use vars qw( @ISA );
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::PredictionTranscript;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );



=head2 fetch_by_dbID

  Arg  1    : int $dbID
              database internal id for a PredictionTranscript
  Function  : Retrieves PredictionTranscript from db with given dbID.
  Returntype: Bio::EnsEMBL::PredictionTranscript
  Exceptions: returns undef when not found
  Caller    : general

=cut



sub fetch_by_dbID {
  my ( $self, $dbID ) = @_;
  my $hashRef;
  my @exons;

  if( !defined $dbID ) {
      $self->throw("Give a prediction_transcript_id");
  }

  my $query = qq {
    SELECT  p.prediction_transcript_id
      , p.contig_id
      , p.contig_start
      , p.contig_end
      , p.contig_strand
      , p.start_phase
      , p.exon_rank
      , p.score
      , p.p_value	
      , p.analysis_id
      , p.exon_count

    FROM prediction_transcript p
    WHERE p.prediction_transcript_id = ?
    ORDER BY p.prediction_transcript_id, p.exon_rank
  };

  my $sth = $self->prepare( $query );
  $sth->execute( $dbID );

  my @res = $self->_ptrans_from_sth( $sth );
  return $res[0];
}



=head2 fetch_by_Slice

  Arg  1    : Bio::EnsEMBL::Slice $slice
              start, end, chromosome and type are used from the slice.
  Function  : retrieves PTs from database which overlap with given slice.
              PTs not fully in the slice are retrieved partially, with exons 
              set to undef. 
  Returntype: listref Bio::EnsEMBL::PredictionTranscript
  Exceptions: none, list can be empty.
  Caller    : Slice or WebSlice

=cut

sub fetch_by_Slice {
   my $self = shift;
   my $slice = shift;

   $self->throw( "Not implemented yet" );
}

sub fetch_by_Slice_and_logic_name {
   my $self = shift;
   my $slice = shift;
   my $logic_name = shift;

   $self->throw( "Not implemented yet" );
}


=head2 fetch_by_Contig

  Arg  1    : Bio::EnsEMBL::RawContig $contig
              Only dbID in Contig is used.
  Function  : returns all PredicitonTranscipts on given contig
  Returntype: listref Bio::EnsEMBL::PredictionTranscript
  Exceptions: none, if there are none, the list is empty.
  Caller    : Bio::EnsEMBL::RawContig->get_genscan_peptides();

=cut


sub fetch_by_Contig {
  my $self = shift;
  my $contig = shift;

  my $query = qq {
    SELECT  p.prediction_transcript_id
      , p.contig_id
      , p.contig_start
      , p.contig_end
      , p.contig_strand
      , p.start_phase
      , p.exon_rank
      , p.score
      , p.p_value	
      , p.analysis_id
      , p.exon_count

    FROM prediction_transcript p
    WHERE p.contig_id = ?
    ORDER BY p.prediction_transcript_id, p.exon_rank
  };

  my $sth = $self->prepare( $query );
  $sth->execute( $contig->dbID );

  my @res = $self->_ptrans_from_sth( $sth );
  return \@res;
}

sub fetch_by_Contig_and_logic_name {
  my $self = shift;
  my $contig = shift;
  my $logic_name = shift;

  my $query = qq {
    SELECT  p.prediction_transcript_id
      , p.contig_id
      , p.contig_start
      , p.contig_end
      , p.contig_strand
      , p.start_phase
      , p.exon_rank
      , p.score
      , p.p_value	
      , p.analysis_id
      , p.exon_count

    FROM prediction_transcript p, analysis a
    WHERE p.contig_id = ? 
    AND a.analysis_id = p.analysis_id
    AND a.logic_name = ?
    ORDER BY p.prediction_transcript_id, p.exon_rank
  };

  my $sth = $self->prepare( $query );
  $sth->execute( $contig->dbID, $logic_name );

  my @res = $self->_ptrans_from_sth( $sth );
  return \@res;
}


=head2 _ptrans_from_sth

  Arg  1    : DBI:st $statement_handle
              an already executed statement handle.
  Function  : Generate PredictionTranscripts from the handle. Obviously 
              this needs to come from a query on prediciton_transcript.
              Needs to be sorted on exon_rank and p.._t.._id.
  Returntype: list Bio::EnsEMBL::PredictionTranscript
  Exceptions: none, list can be empty
  Caller    : internal

=cut


sub _ptrans_from_sth {
  my $self = shift;
  my $sth = shift;

  my $analysis;
  my $pre_trans = undef; 
  my $pre_trans_id = undef;
  my @result = ();

  while( my $hashRef = $sth->fetchrow_hashref() ) {
    if(( ! defined $pre_trans  ) ||
       ( $pre_trans_id != $hashRef->{'prediction_transcript_id'} )) {

      $pre_trans = Bio::EnsEMBL::PredictionTranscript->new(); 
      $pre_trans_id = $hashRef->{'prediction_transcript_id'};
      $pre_trans->dbID( $pre_trans_id );
      $pre_trans->adaptor( $self );
      my $anaAdaptor = $self->db()->get_AnalysisAdaptor();
      $analysis = $anaAdaptor->fetch_by_dbID( $hashRef->{'analysis_id'} );
      $pre_trans->analysis( $analysis );
      $pre_trans->set_exon_count( $hashRef->{'exon_count'} );
      push( @result, $pre_trans );
    }

    my $exon = $self->_new_Exon_from_hashRef( $hashRef );
    $pre_trans->add_Exon( $exon, $hashRef->{'exon_rank'} );
  }
  return @result;
}


=head2 _new_Exon_from_hashRef

  Arg  1    : hashref $exon_attributes
              Data from a line in prediction_transcript
  Function  : Creates an Exon from the data
  Returntype: Bio::EnsEMBL::Exon
  Exceptions: none
  Caller    : internal

=cut



sub _new_Exon_from_hashRef {
  my $self = shift;
  my $hashRef = shift;
  
  my $exon = Bio::EnsEMBL::Exon->new();
  my $contig_adaptor = $self->db()->get_RawContigAdaptor();
  
  my $contig = Bio::EnsEMBL::RawContig->new
    ( $hashRef->{'contig_id'}, $contig_adaptor );
  
  $exon->start( $hashRef->{'contig_start'} );
  $exon->end( $hashRef->{'contig_end'} );
  $exon->strand( $hashRef->{'contig_strand'} );
  $exon->phase( $hashRef->{start_phase} );
  
  $exon->contig( $contig );
  $exon->attach_seq( $contig->seq() );
  $exon->ori_start( $exon->start );
  $exon->ori_end( $exon->end );
  $exon->ori_strand( $exon->strand );
  
  # does exon not have score?
  $exon->score( $hashRef->{'score'} );
  $exon->p_value( $hashRef->{'p_value'} );
  
  return $exon;
}





=head2 store

  Arg  1    : Bio::EnsEMBL::PredictionTranscript $pt
  Function  : Stores given $pt in database. Puts dbID and Adaptor into $pt 
              object. Returns the dbID.
  Returntype: int
  Exceptions: on wrong argument type
  Caller    : general

=cut

sub store {
  my ( $self, $pre_trans ) = @_;

  if( ! $pre_trans->isa('Bio::EnsEMBL::PredictionTranscript') ) {
    $self->throw("$pre_trans is not a EnsEMBL PredictionTranscript - not dumping!");
  }

  if( $pre_trans->dbID && $pre_trans->adaptor == $self ) {
    $self->warn("Already stored");
  }


  my $exon_sql = q{
       INSERT into prediction_transcript ( prediction_transcript_id, exon_rank, contig_id, 
                                           contig_start, contig_end, contig_strand, 
                                           start_phase, score, p_value,
                                           analysis_id, exon_count )
		 VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )
	       };

  my $exonst = $self->prepare($exon_sql);

  my $exonId = undef;

  my @exons = $pre_trans->get_all_Exons();
  my $dbID = undef;
  my $rank = 1;
  
  for my $exon ( @exons ) {
    if( ! defined $exon ) { $rank++; next; }
    
    my $contig_id = $exon->contig->dbID();
    my $contig_start = $exon->start();
    my $contig_end = $exon->end();
    my $contig_strand = $exon->strand();
    
    my $start_phase = $exon->phase();
    my $end_phase = $exon->end_phase();

    # this is only in PredictionExon
    my $score = $exon->score();
    my $p_value = $exon->p_value();
    #print "storing exon with pvalue ".$exon->p_value." and score ".$exon->score."\n";
    my $analysis = $pre_trans->analysis->dbID;

    if( $rank == 1 ) {
      $exonst->execute( undef, 1, $contig_id, $contig_start, $contig_end, $contig_strand,
			$start_phase, $score, $p_value, $analysis, scalar( @exons) );
      $dbID =   $exonst->{'mysql_insertid'};
    } else {
      $exonst->execute( $dbID, $rank, $contig_id, $contig_start, $contig_end, $contig_strand,
			$start_phase, $score, $p_value, $analysis, scalar( @exons ) );
    }
    $rank++;
  }

  $pre_trans->dbID( $dbID );
  $pre_trans->adaptor( $self );
  
  return $dbID;
}



=head2 remove

  Arg  1    : Bio::EnsEMBL::PredictionTranscript $pt
  Function  : removes given $pt from database. Expects access to
              internal db via $pt->{'dbID'} to set it undef.
  Returntype: none
  Exceptions: none
  Caller    : general

=cut


sub remove {
  my $self = shift;
  my $pre_trans = shift;
  
  if ( ! defined $pre_trans->dbID() ) {
    return;
  }

  my $sth = $self->prepare( "delete from prediction_transcript where prediction_transcript_id = ?" );
  $sth->execute( $pre_trans->dbID );

  # uhh, didnt know another way of resetting to undef ...
  $pre_trans->{dbID} = undef;
  $pre_trans->{adaptor} = undef;
}







1;

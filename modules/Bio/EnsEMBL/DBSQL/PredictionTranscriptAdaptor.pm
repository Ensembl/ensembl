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
use Bio::EnsEMBL::PredictionExon;
use Bio::EnsEMBL::PredictionTranscript;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );


=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   : 
 Function: fetches a PredictionTranscript by its dbID
 Example : 
 Returns : 
 Args    : 

=cut

# returns list of exons or maybe empty list (gene not known)

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
      , p.end_phase
      , p.rank
      , p.score
      , p.analysis_id
      , c.id cid
    FROM prediction_transcript p
      , contig c
    WHERE p.prediction_transcript_id = ?
      AND p.contig_id = c.internal_id
    ORDER BY p.prediction_transcript_id, p.rank
  };

  my $sth = $self->prepare( $query );
  $sth->execute();
  my $analysis = undef;

  while( $hashRef = $sth->fetchrow_hashref() ) {
    if( ! defined $analysis ) {
      $analysis = $db->get_AnalysisAdaptor->fetch_by_dbID( $hashRef->{'analysis_id'} );
    }
    my $exon = $self->_exon_from_sth( $sth, $hashRef );
    push( @exons, $exon );
  }

  delete $self->{rchash};
  
  my $pre_trans = Bio::EnsEMBL::PredictionTranscript->new();
  $pre_trans->analysis( $analysis );
  # add the exons somehow
  my $transl = Bio::EnsEMBL::Translation->new( );
  $transl->start_exon( $exons[0] );
  $transl->start( 1 );
  $transl->end_exon( $exons[$#exons] );
  $transl->end( $exons[$#exons]->length() );

  $pre_trans->translation( $transl );

  return;
}


sub fetch_by_assembly_location {
}

sub fetch_by_contig {
}


sub _new_Exon_from_hashRef {
   my $self = shift;
   my $hashRef = shift;

   my $exon = Bio::EnsEMBL::Exon->new();

   #  if exon doesnt support score we need PredictionExons which should
   # extend Exon
   #   my $exon = Bio::EnsEMBL::PredictionExon->new();
   $exon->start( $hashRef->{'contig_start'} );
   $exon->end( $hashRef->{'contig_end'} );
   $exon->strand( $hashRef->{'contig_strand'} );
   $exon->phase( $hashRef->{start_phase} );
   $exon->end_phase( $hashRef->{end_phase} );
   $exon->adaptor($self);
   
   if( !exists $self->{rchash}{$hashRef->{'contig_id'}} ) {
     $self->{rchash}{$hashRef->{contig_id}} = $self->db->get_Contig($hashRef->{'cid'});
   }

   $exon->attach_seq($self->{rchash}{$hashRef->{'contig_id'}}->primary_seq);
   $exon->contig( $self->{rchash}{$hashRef->{'contig_id'}} );
   $exon->seqname($hashRef->{'cid'});
   $exon->ori_start( $exon->start );
   $exon->ori_end( $exon->end );
   $exon->ori_strand( $exon->strand );

   # does exon not have score?
   $exon->score( $hashRef->{'score'} );

# consider the following ...
#   $exon->dbID( $hashRef->{'rank'} )
   
  return $exon;
}





=head2 store

 Title   : store
 Usage   : $ptransAdaptor->store($predicitionTranscriptObject)
 Function: Stores the thing
 Example : 
 Returns : 
 Args    : predictionTranscript object

=cut

sub store {
  my ( $self, $pre_trans ) = @_;

  if( ! $pre_trans->isa('Bio::EnsEMBL::Prediction_Transcript') ) {
    $self->throw("$pre_trans is not a EnsEMBL PredictionTranscript - not dumping!");
  }

  if( $pre_trans->dbID && $pre_trans->adaptor == $self ) {
    $self->warn("Already stored");
  }


  my $exon_sql = q{
       INSERT into prediction_transcript ( prediction_transcript_id, exon_rank, contig_id, contig_start, contig_end, contig_strand, start_phase, end_phase, score, analysis_id )
		 VALUES ( ?, ?, ?, ?, ?, ?, ?,?, ? )
	       };
  my $exonst = $self->prepare($exon_sql);

  my $exonId = undef;

  my @exons = $pre_trans->get_all_exons();
  my $dbID = undef;
  my $rank = 1;

  for my $exon ( @exons ) {

    my $contig_id = $exon->contig->dbID();
    my $contig_start = $exon->start();
    my $contig_end = $exon->end();
    my $contig_strand = $exon->strand();
    
    my $start_phase = $exon->phase();
    my $end_phase = $exon->end_phase();

    # this is only in PredictionExon
    my $score = $exon->score();

    my $analysis = $pre_trans->analysis->dbID;

    if( $rank = 1 ) {
      $exonst->execute( undef, 1, $contig_id, $contig_start, $contig_end, $contig_strand,
			$start_phase, $end_phase, $score, $analysis );
      $dbID =   $exonst->{'mysql_insertid'};
    } else {
      $exonst->execute( $dbID, $rank, $contig_id, $contig_start, $contig_end, $contig_strand,
			$start_phase, $end_phase, $score, $analysis );
    }
    $rank++;
  }

  $pre_trans->dbID( $dbID );
  $pre_trans->adaptor( $self );
  
  return $dbID;
}


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
}







1;

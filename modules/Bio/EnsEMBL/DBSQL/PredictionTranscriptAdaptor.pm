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

#get a prediction transcript adaptor from the database
$pta = $database_adaptor->get_PredictionTranscriptAdaptor();

#get a slice on a region of chromosome 1
$sa = $database_adaptor->get_SliceAdaptor();
$slice = $sa->fetch_by_chr_start_end('1', 100000, 200000);

#get all the prediction transcripts from the slice region
$prediction_transcripts = @{$pta->fetch_all_by_Slice($slice)};

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

  Arg [1]    : int $dbID
               database internal id for a PredictionTranscript 
  Example    : $pt = $prediction_transcript_adaptor->fetch_by_dbID($dbID);
  Description: Retrieves PredictionTranscript from db with given dbID. 
  Returntype : Bio::EnsEMBL::PredictionTranscript in contig coords
  Exceptions : none
  Caller     : general

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

  return $self->_ptrans_from_sth( $sth )->[0];
}


=head2 fetch_all_by_RawContig

  Arg [1]    : Bio::EnsEMBL::RawContig $contig
               The contig to retrieve prediction transcripts from
  Arg [2]    : (optional) string $logic_name
               the type of analysis performed on objects which should be
               retrieved
  Example    : @pts = @{$pt_adptr->fetch_all_by_RawContig($contig,'Genscan')};
  Description: Retrieves prediction transcripts from a contig
  Returntype : listref of Bio::EnsEMBL::PredictionTranscripts in contig coords
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_RawContig {
  my ($self, $contig, $logic_name) = @_;

  my $constraint = undef;
  
  if($logic_name){
    my $analysis  = 
      $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);

    unless($analysis && $analysis->dbID()) {
      $self->warn("no analysis with logic_name '$logic_name' exists");
      return ();
    }

    $constraint = " analysis_id = ".$analysis->dbID;
  }

  return $self->fetch_all_by_RawContig_constraint($contig, $constraint);
}


=head2 fetch_all_by_RawContig_constraint

  Arg [1]    : Bio::EnsEMBL::RawContig $contig
               The contig to obtain prediction transcripts from
  Arg [2]    : (optional) string $constraint
               the limiting SQL to form the where clause of the the database
               query
  Example    : $pts = $pta->fetch_all_by_RawContig_constraint($contig, 
                                                            'analysis_id = 2');
  Description: returns all PredicitonTranscipts on given contig 
  Returntype : listref of Bio::EnsEMBL::PredictionTranscript in contig coords 
  Exceptions : none, if there are none, the list is empty.
  Caller     : ?

=cut

sub fetch_all_by_RawContig_constraint {
  my $self = shift;
  my $contig = shift;
  my $constraint = shift;

  unless(defined $contig && $contig->isa('Bio::EnsEMBL::RawContig')){
    $self->throw("contig arg must be a Bio::EnsEMBL::RawContig");
    return undef;
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
    WHERE p.contig_id = ?
   };

  if($constraint){
    $query .= " and ".$constraint;
  }

  $query .= " order by p.prediction_transcript_id, p.exon_rank";
  #print $query."\n";
  my $sth = $self->prepare( $query );
  $sth->execute( $contig->dbID() );

  return $self->_ptrans_from_sth( $sth );
}


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice of the region to obtain prediction transcripts from
  Arg [2]    : (optional) string $logic_name
               the type of analysis performed on the prediction transcripts to 
               obtain
  Example    : @pts = $pta->fetch_all_by_Slice($slice,'Genscan');
  Description: returns all PredicitonTranscipts on the region of the slice 
  Returntype : listref of Bio::EnsEMBL::PredictionTranscript in slice coords 
  Exceptions : none, if there are none, the list is empty.
  Caller     : ?

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $logic_name) = @_;

  my $constraint = undef;

  if($logic_name){
    my $analysis  = 
      $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);

    unless($analysis->dbID()) {
      $self->warn("PredictionTranscriptAdaptor::fetch_all_by_slice: " .
		  "no analysis with logic name $logic_name exists\n");
      return ();
    }
    $constraint = " analysis_id = ".$analysis->dbID;
  }
  
  my $results = 
    $self->fetch_all_by_assembly_location_constraint($slice->chr_start, 
						 $slice->chr_end, 
						 $slice->chr_name, 
						 $slice->assembly_type, 
						 $constraint);

  my @out = ();

 GENE: foreach my $transcript(@$results) {
    my $exon_count = 1;
    my $pred_t = Bio::EnsEMBL::PredictionTranscript->new();
    $pred_t->dbID($transcript->dbID);
    $pred_t->adaptor($self);
    $pred_t->analysis($transcript->analysis);
    $pred_t->set_exon_count($transcript->get_exon_count);
    my $exons = $transcript->get_all_Exons;
    my @sorted_exons;
    if($exons->[0]->strand == 1){
      @sorted_exons = sort{$a->start <=> $b->start} @$exons;
    }else{
      @sorted_exons = sort{$b->start <=> $a->start} @$exons;
    }
    my $contig = $sorted_exons[0]->contig;
  EXON:foreach my $e(@sorted_exons){ 
      my ( $start, $end, $exon );

      if( $slice->strand == 1 ) {
	$start = ($e->start - ($slice->chr_start - 1));
	$end = ($e->end - ($slice->chr_start - 1));
	$exon = $self->_new_Exon($start, $end, $e->strand, 
				    $e->phase, $e->score, $e->p_value, $contig);
      } else {
	$start = $slice->chr_end() - $e->end() + 1;
	$end = $slice->chr_end() - $e->start() + 1;
	$exon = $self->_new_Exon($start, $end, -$e->strand, 
				    $e->phase, $e->score, $e->p_value, $contig);
      }      
      $pred_t->add_Exon( $exon, $exon_count );
      $exon_count++;
      push(@out, $pred_t);
    }
  } 
  return \@out;
}


=head2 fetch_all_by_assembly_location

  Arg [1]    : int $chr_start 
               the start of the region to obtain prediction transcripts from
               in chromosomal coordinates
  Arg [2]    : int $chr_end
               the end of the region to obtain prediction transcripts from
               in chromosomal coordinates
  Arg [3]    : string $chr
               the chromosome to obtain prediction transcripts from
  Arg [4]    : string $type
               the type of assembly to use to obtain transcripts from
  Arg [5]    : (optional) string $logic_name
               the type of analysis used to construct the prediction 
               transcripts
  Example    : @pts = $pta->fetch_all_by_assembly_location(1,20000, '2', 
							  'NCBI30', 'Genscan');
  Description: Retrieves prediction transcripts from a region of the assembly 
  Returntype : listref of Bio::EnsEMBL::PredictionTranscripts in chromo coords
  Exceptions : none
  Caller     : fetch_all_by_Slice

=cut

sub fetch_all_by_assembly_location{
  my ($self, $chr_start, $chr_end, $chr, $type, $logic_name) = @_;

  my $constraint = undef;

  if($logic_name){
    my $analysis = 
      $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
    $constraint = " analysis_id = ".$analysis->dbID;

    unless($analysis && $analysis->dbID()) {
      $self->warn("no analysis exists for logic_name $logic_name");
      return ();
    }
  }
  
  return $self->fetch_all_by_assembly_location_constraint($chr_start, 
						      $chr_end, 
						      $chr, 
						      $type, 
						      $constraint);

}


=head2 fetch_all_by_assembly_location_constraint

  Arg [1]    : int $chr_start 
               the start of the region to obtain prediction transcripts from
               in chromosomal coordinates
  Arg [2]    : int $chr_end
               the end of the region to obtain prediction transcripts from
               in chromosomal coordinates
  Arg [3]    : string $chr
               the chromosome to obtain prediction transcripts from
  Arg [4]    : string $type
               the type of assembly to use to obtain transcripts from
  Arg [5]    : (optional) string $contraint
               SQL constraint to limit the values returned (inserted in 
               WHERE clause).
  Example    : $pts = $pta->fetch_all_by_assembly_locationconstraint(1,20000, 
						'2', 'NCBI30','analysis_id=2');
  Description: Retrieves prediction transcripts from a region of the assembly 
  Returntype : listref of Bio::EnsEMBL::PredictionTranscripts in chromo coords
  Exceptions : thrown if $chr_start or $chr_end are not numbers
  Caller     : fetch_all_by_assembly_location

=cut

sub fetch_all_by_assembly_location_constraint {
  my ($self, $chr_start, $chr_end, $chr, $type, $constraint) = @_;

  if( !defined $type ) {
    $self->throw("Assembly location must be start,end,chr,type");
  }
  
  if( $chr_start !~ /^\d/ || $chr_end !~ /^\d/ ) {
    $self->throw("start/end must be numbers not $chr_start,$chr_end " .
		 "(have you typed the location in the right way around" .
		 " - start,end,chromosome,type)?");
  }
  
  my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type($type);
  
  my @cids = $mapper->list_contig_ids($chr, $chr_start ,$chr_end);
  my %ana;
  my $cid_list = join(',',@cids);
  
  my $sql = qq {
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
    WHERE
   };

  $sql .= "contig_id in($cid_list) ";

  if($constraint){
    $sql .= " and $constraint";
  }

  my $sth = $self->prepare($sql);
  $sth->execute;

  my $results = $self->_ptrans_from_sth($sth);
  my @out = ();
  GENE: foreach my $transcript(@$results){
      my $exon_count = 1;
      my $pred_t = Bio::EnsEMBL::PredictionTranscript->new();
      $pred_t->dbID($transcript->dbID);
      $pred_t->adaptor($self);
      $pred_t->analysis($transcript->analysis);
      $pred_t->set_exon_count($transcript->get_exon_count);
      my $exons = $transcript->get_all_Exons;
      my @sorted_exons;
      if($exons->[0]->strand == 1){
	@sorted_exons = sort{$a->start <=> $b->start} @$exons;
      }else{
	@sorted_exons = sort{$b->start <=> $a->start} @$exons;
      }
      my $contig = $sorted_exons[0]->contig;
    EXON:foreach my $e(@sorted_exons){
	my @coord_list = 
	  $mapper->map_coordinates_to_assembly($e->contig->dbID, $e->start, 
					       $e->end, $e->strand, 
					       "rawcontig");
	if( scalar(@coord_list) > 1 ) {
	  #$self->warn("maps to ".scalar(@coord_list)." .
	  #             "coordinate objs not all of feature will be on " .
	  #             "golden path skipping\n");
	  next GENE;
	}
	
	if($coord_list[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
	  #$self->warn("this feature is on a part of $contig_id which isn't " .
	  #             "on the golden path skipping");
	  next GENE;
	}
	if(!($coord_list[0]->start >= $chr_start) ||
	   !($coord_list[0]->end <= $chr_end)) {
	  next GENE;
	}
	my $exon = $self->_new_Exon($coord_list[0]->start, $coord_list[0]->end,
				    $coord_list[0]->strand, $e->phase, 
				    $e->score, $e->p_value, $contig);
	$pred_t->add_Exon( $exon, $exon_count );
	$exon_count++;
      }
      push(@out, $pred_t);
    }

  return \@out;
}


=head2 _ptrans_from_sth

  Arg [1]    : DBI:st $statement_handle
               an already executed statement handle. 
  Example    : none
  Description: PRIVATE 
               Generate PredictionTranscripts from the handle. Obviously 
               this needs to come from a query on prediciton_transcript.
               Needs to be sorted on exon_rank and p.._t.._id. 
  Returntype : list Bio::EnsEMBL::PredictionTranscript in contig coords 
  Exceptions : none
  Caller     : internal

=cut

sub _ptrans_from_sth {
  my $self = shift;
  my $sth = shift;

  my $analysis;
  my $pre_trans = undef; 
  my $pre_trans_id = undef;
  my @result = ();
  my $count = 0;
  my $exon_count = 0;
  while( my $hashRef = $sth->fetchrow_hashref() ) {
    if(( ! defined $pre_trans  ) ||
       ( $pre_trans_id != $hashRef->{'prediction_transcript_id'} )) {
      $count++;
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
    $exon_count++;
  }
  #print "have created ".$count." transcripts and ".$exon_count." exons\n";
  return \@result;
}


=head2 _new_Exon_from_hashRef

  Arg [1]    : hashref $exon_attributes
               Data from a line in prediction_transcript 
  Example    : none
  Description: PRIVATE Creates an Exon from the data 
  Returntype : Bio::EnsEMBL::Exon in contig coords
  Exceptions : none
  Caller     : internal

=cut

sub _new_Exon_from_hashRef {
  my $self = shift;
  my $hashRef = shift;
  
  my $exon = Bio::EnsEMBL::Exon->new();
  my $contig_adaptor = $self->db()->get_RawContigAdaptor();
  
  my $contig = $contig_adaptor->fetch_by_dbID($hashRef->{'contig_id'});
  
  $exon->start( $hashRef->{'contig_start'} );
  $exon->end( $hashRef->{'contig_end'} );
  $exon->strand( $hashRef->{'contig_strand'} );
  $exon->phase( $hashRef->{start_phase} );
  $exon->end_phase( ($exon->end - $exon->start + 1 + $exon->phase) % 3 );
  
  $exon->contig( $contig );
  $exon->attach_seq( $contig );
  $exon->ori_start( $exon->start );
  $exon->ori_end( $exon->end );
  $exon->ori_strand( $exon->strand );
  
  # does exon not have score?
  $exon->score( $hashRef->{'score'} );
  $exon->p_value( $hashRef->{'p_value'} );
  
  return $exon;
}



=head2 _new_Exon

  Arg [1]    : int $start
  Arg [2]    : int $end
  Arg [3]    : int $strand
  Arg [4]    : int $phase
  Arg [5]    : float $score
  Arg [6]    : float $pvalue
  Arg [7]    : Bio::EnsEMBL::Contig $contig
  Example    : none 
  Description: PRIVATE creates a new exon from args and returns it
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : internal

=cut

sub _new_Exon{
  my ($self, $start, $end, $strand, $phase, $score, $pvalue, $contig) = @_; 
  my $exon = Bio::EnsEMBL::Exon->new();
  
  $exon->start( $start);
  $exon->end( $end );
  $exon->strand( $strand );
  $exon->phase( $phase );
  
  $exon->contig( $contig );
  $exon->attach_seq( $contig );
  $exon->ori_start( $start );
  $exon->ori_end( $end );
  $exon->ori_strand( $strand );
  
  # does exon not have score?
  $exon->score( $score );
  $exon->p_value( $pvalue );
  
  return $exon;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::PredictionTranscript $pre_trans 
  Example    : $prediction_transcript_adaptor->store($pre_trans);
  Description: Stores given $pt in database. Puts dbID and Adaptor into $pt 
               object. Returns the dbID. 
  Returntype : int 
  Exceptions : on wrong argument type 
  Caller     : general 

=cut

sub store {
  my ( $self, $pre_trans ) = @_;

  if( ! $pre_trans->isa('Bio::EnsEMBL::PredictionTranscript') ) {
    $self->throw("$pre_trans is not a EnsEMBL PredictionTranscript " 
		 . "- not dumping!");
  }

  if( $pre_trans->dbID && $pre_trans->adaptor == $self ) {
    $self->warn("Already stored");
  }


  my $exon_sql = q{
    INSERT INTO prediction_transcript ( prediction_transcript_id, exon_rank, 
					contig_id, contig_start, contig_end, 
					contig_strand, start_phase, score, 
					p_value, analysis_id, exon_count )
    VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )
  };

  my $exonst = $self->prepare($exon_sql);

  my $exonId = undef;

  my $exons = $pre_trans->get_all_Exons();
  my $dbID = undef;
  my $rank = 1;
  
  for my $exon ( @$exons ) {
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

    my $analysis = $pre_trans->analysis->dbID;

    if( $rank == 1 ) {
      $exonst->execute( undef, 1, $contig_id, $contig_start, 
			$contig_end, $contig_strand,
			$start_phase, $score, $p_value, $analysis, 
			scalar( @{$exons} ));
      $dbID = $exonst->{'mysql_insertid'};
    } else {
      $exonst->execute( $dbID, $rank, $contig_id, $contig_start, 
			$contig_end, $contig_strand,
			$start_phase, $score, $p_value, $analysis, 
			scalar( @{$exons} ) );
    }
    $rank++;
  }

  $pre_trans->dbID( $dbID );
  $pre_trans->adaptor( $self );
  
  return $dbID;
}



=head2 remove

  Arg [1]    : Bio::EnsEMBL::PredictionTranscript $pt 
  Example    : $prediction_transcript_adaptor->remove($pt);
  Description: removes given prediction transcript $pt from database. 
  Returntype : none
  Exceptions : none 
  Caller     : general

=cut

sub remove {
  my $self = shift;
  my $pre_trans = shift;
  
  if ( ! defined $pre_trans->dbID() ) {
    return;
  }

  my $sth = $self->prepare( "DELETE FROM prediction_transcript 
                             WHERE prediction_transcript_id = ?" );
  $sth->execute( $pre_trans->dbID );

  # uhh, didnt know another way of resetting to undef ...
  $pre_trans->{dbID} = undef;
  $pre_trans->{adaptor} = undef;
}



=head2 fetch_by_Contig

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_RawContig instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Contig {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Contig has been renamed fetch_all_by_RawContig\n" . caller);

  return $self->fetch_all_by_RawContig(@args);
}


=head2 fetch_by_Contig_constraint

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_by_RawContig_constraint instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Contig_constraint {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Contig_constraint has been renamed fetch_by_RawContig_constraint\n" . caller);

  return $self->fetch_by_RawContig_constraint(@args);
}


=head2 fetch_by_Slice

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_Slice instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Slice {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Slice has been renamed fetch_all_by_Slice\n" . caller);

  return $self->fetch_all_by_Slice(@args);
}


=head2 fetch_by_assembly_location

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_assembly_location instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_assembly_location {
  my ($self, @args) = @_;

  $self->warn("fetch_by_assembly_location has been renamed fetch_all_by_assembly_location\n" . caller);

  return $self->fetch_all_by_assembly_location(@args);
}


=head2 fetch_by_assembly_location_constraint

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_assembly_location_constraint instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_assembly_location_constraint {
  my ($self, @args) = @_;

  $self->warn("fetch_by_assembly_location_constraint has been renamed fetch_all_by_assembly_location_constraint\n" . caller);

  return $self->fetch_all_by_assembly_location_constraint(@args);
}



1;

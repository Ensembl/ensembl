# EnsEMBL Exon reading writing adaptor for mySQL
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

Email questions to the EnsEMBL developer list: <ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor;

use vars qw( @ISA );
use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::PredictionTranscript;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor );


=head2 _tables

  Arg [1]    : none
  Example    : none
  Description: Implements abstract superclass method to define the table used
               to retrieve prediction transcripts from the database
  Returntype : string
  Exceptions : none
  Caller     : generic_fetch

=cut

sub _tables {
  my $self = shift;

  return ['prediction_transcript', 'p'];
}



=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: Implements abstract superclass method to define the columns
               retrieved in database queries used to create prediction 
               transcripts.
  Returntype : list of strings
  Exceptions : none
  Caller     : generic_fetch

=cut

sub _columns {
  my $self = shift;

  return qw( p.prediction_transcript_id
       p.contig_id
       p.contig_start
       p.contig_end
       p.contig_strand
       p.start_phase
       p.exon_rank
       p.score
       p.p_value	
       p.analysis_id
       p.exon_count);
}



=head2 _final_clause

  Arg [1]    : none
  Example    : none
  Description: Overrides superclass method to provide an additional table
               joining coinstraint before the SQL query is performed.
  Returntype : string
  Exceptions : none
  Caller     : generic_fetch

=cut

sub _final_clause {
  my $self = shift;
 
  return  'order by p.prediction_transcript_id, p.exon_rank';
}

=head2 fetch_by_stable_id

Arg [1]    : string $stable_id
             The stable id of the transcript to retrieve
Example    : $trans = $trans_adptr->fetch_by_stable_id('3.10.190');
Description: Retrieves a prediction transcript via its stable id.  Note that
             the stable id is not actually stored in the database and is 
             calculated upon retrieval as the contig.start.end of the 
             prediction transcript
Returntype : Bio::EnsEMBL::PredictionTranscript
Caller     : general

=cut

sub fetch_by_stable_id {
  my $self = shift;
  my $id = shift;

  $self->throw('ID argument is required') if(!$id);

  #strip the start/end off the end of the id to get the contig name
  my @components = split(/\./, $id);
  pop(@components);
  pop(@components);

  my $contig_name = join('.', @components);
  my $contig = $self->db->get_RawContigAdaptor->fetch_by_name($contig_name);

  $self->throw("unknown contig [$contig_name] in identifier [$id]") 
    if(!$contig);
  
  #
  # Retrieve all the prediction transcripts on the contig it is known to be on.
  # Return the one found with the matching 'id'
  #
  foreach my $pt (@{$self->fetch_all_by_RawContig($contig)}) {
    return $pt if($id eq $pt->stable_id);
  }

  #prediction transcript does not exist
  $self->throw("No prediction transcipt exists with stable id [$id]");
}


=head2 _objs_from_sth

  Arg [1]    : DBI:st $sth 
               An executed DBI statement handle
  Arg [2]    : (optional) Bio::EnsEMBL::Mapper $mapper 
               An mapper to be used to convert contig coordinates
               to assembly coordinates.
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice to map the prediction transcript to.   
  Example    : $p_transcripts = $self->_objs_from_sth($sth);
  Description: Creates a list of Prediction transcripts from an executed DBI
               statement handle.  The columns retrieved via the statement 
               handle must be in the same order as the columns defined by the
               _columns method.  If the slice argument is provided then the
               the prediction transcripts will be in returned in the coordinate
               system of the $slice argument.  Otherwise the prediction 
               transcripts will be returned in the RawContig coordinate system.
  Returntype : reference to a list of Bio::EnsEMBL::PredictionTranscripts
  Exceptions : none
  Caller     : superclass generic_fetch

=cut

sub _objs_from_sth {
  my ($self, $sth, $mapper, $slice) = @_;
  
  my @out = ();
  
  my ($prediction_transcript_id, 
      $contig_id, $contig_start, $contig_end, $contig_strand,
      $start_phase, $exon_rank, $score, $p_value, $analysis_id,
      $exon_count );

  $sth->bind_columns(\$prediction_transcript_id, 
		    \$contig_id, \$contig_start, \$contig_end, \$contig_strand,
		    \$start_phase, \$exon_rank, \$score, \$p_value, 
		    \$analysis_id,\$exon_count);

  my $rca = $self->db->get_RawContigAdaptor;
  my $aa  = $self->db->get_AnalysisAdaptor;
  
  my ($analysis, $contig, $pre_trans, $ptid, $on_slice_flag, $last_end,
      $chr, $start, $end, $strand, 
      $slice_start, $slice_end, $slice_strand,
      $exon, $exon_start, $exon_end, $exon_strand,
      $stable_start, $stable_end, $stable_ctg,
      $transcript_slice_start, $transcript_slice_end );
  my (%analysis_hash, %contig_hash);

  if($slice) {
    $slice_start  = $slice->chr_start;
    $slice_end    = $slice->chr_end;
    $slice_strand = $slice->strand;
  }

  $on_slice_flag = 0;

  my $prev_exon;
  my $already_merged;
  
  while($sth->fetch) {
    #create a new transcript for each new prediction transcript id
    unless(defined $pre_trans && $ptid == $prediction_transcript_id) {
      $pre_trans = Bio::EnsEMBL::PredictionTranscript->new;

      $ptid = $prediction_transcript_id;
      $pre_trans->dbID($ptid);
      
      unless($analysis = $analysis_hash{$analysis_id}) {
	$analysis = $aa->fetch_by_dbID($analysis_id);
	$analysis_hash{$analysis_id} = $analysis;
      }
      
      $pre_trans->analysis($analysis);
      $pre_trans->set_exon_count($exon_count);
  
      $prev_exon = undef;
      $already_merged = 0;

      if(@out) {
	#throw away last pt if no exons or introns were on the slice
	if($slice && ( $transcript_slice_end < 1 || 
		       $transcript_slice_start > $slice->length() )) {
	  pop @out;
	} else {
	  #set the stable_id of the previous prediction
	  $out[$#out]->stable_id("$stable_ctg.$stable_start.$stable_end");
	}
      }
      
      push( @out, $pre_trans );

      #reset values used for last predtrans
      $stable_start = -1;
      $stable_end   = -1;
      $stable_ctg = '';

      $transcript_slice_end = undef;
      $transcript_slice_start = undef;
    }

    #recalculate stable id values
    if($stable_start == -1 || $contig_start < $stable_start) {
      $stable_start = $contig_start;
    }
    if($contig_end > $stable_end) {
      $stable_end = $contig_end;
    }
    unless($contig = $contig_hash{$contig_id}) {
      $contig = $rca->fetch_by_dbID($contig_id);
      $contig_hash{$contig_id} = $contig;
    }
    $stable_ctg = $contig->name;

    if($slice) {
      #a slice was passed in so we want slice coords

      #convert contig coords to assembly coords
      ($chr, $start, $end, $strand) = 
	$mapper->fast_to_assembly($contig_id, $contig_start,
				  $contig_end, $contig_strand);
      
      #if mapped to gap skip
      next unless(defined $start);

      
      #convert to slice coordinates
      if($slice_strand == -1) {
	$exon_start  = $slice_end - $end   + 1;
	$exon_end    = $slice_end - $start + 1;
	$exon_strand = $strand * -1;

	#merge adjacent exons into a single exon
	if($prev_exon && $prev_exon->start == $exon->end + 1) {
	  $exon->end($prev_exon->end);
	  $already_merged++;
	}

      } else {
	$exon_start  = $start - $slice_start + 1;
	$exon_end    = $end   - $slice_start   + 1;
	$exon_strand = $strand;

	#merge adjacent exons into a single exon
	if($prev_exon && $exon->start == $prev_exon->end +1) {
	  $exon->start($prev_exon->start);
	  $already_merged++;
	}
      }   

      if( !defined $transcript_slice_start || 
	  $transcript_slice_start > $exon_start ) {
	$transcript_slice_start = $exon_start;
      }
      
      if( ! defined $transcript_slice_end ||
	  $transcript_slice_end < $exon_end ) {
	$transcript_slice_end = $exon_end;
      }
      #use slice as the contig instead of the raw contig
      $contig = $slice;
    } else {
      #we just want plain old contig coords
      $exon_start =  $contig_start;
      $exon_end   =  $contig_end;
      $exon_strand = $contig_strand;
    }

    #create an exon and add it to the prediction transcript
    $exon = Bio::EnsEMBL::Exon->new_fast($contig, 
					 $exon_start, 
					 $exon_end,
					 $exon_strand);
    $exon->phase( $start_phase );
    $exon->end_phase( ($exon_end - $exon_start + 1 + $start_phase) % 3 );
    $exon->score( $score );
    $exon->p_value( $p_value );

    $prev_exon = $exon;
    $pre_trans->add_Exon($exon, $exon_rank - $already_merged);
  }
  
  #throw away last  pred_transcript if it had no exons overlapping the slice
  if(@out) {
    if($slice && ( $transcript_slice_end < 1 || 
		   $transcript_slice_start > $slice->length() )) {
      pop @out;
    } else {
      #set the stable id of the last prediction transcript
      $out[$#out]->stable_id("$stable_ctg.$stable_start.$stable_end");
    }
  }

  return \@out;
}



=head2 store

  Arg [1]    : list of Bio::EnsEMBL::PredictionTranscript @pre_transcripts 
  Example    : $prediction_transcript_adaptor->store(@pre_transcripts);
  Description: Stores a list of given prediction transcripts in database. 
               Puts dbID and Adaptor into each object stored object.
  Returntype : none
  Exceptions : on wrong argument type 
  Caller     : general 

=cut

sub store {
  my ( $self, @pre_transcripts ) = @_;

  my $exon_sql = q{
      INSERT INTO prediction_transcript ( prediction_transcript_id, exon_rank, 
					  contig_id, contig_start, contig_end, 
					  contig_strand, start_phase, score, 
					  p_value, analysis_id, exon_count )
	VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )
      };

  my $exonst = $self->prepare($exon_sql);

  foreach my $pre_trans (@pre_transcripts) {
    if( ! $pre_trans->isa('Bio::EnsEMBL::PredictionTranscript') ) {
      $self->throw("$pre_trans is not a EnsEMBL PredictionTranscript " 
		   . "- not dumping!");
    }
    
    if( $pre_trans->dbID && $pre_trans->adaptor == $self ) {
      $self->warn("Already stored");
    }
        
    my $exonId = undef;    
    my @pt_exons = @{$pre_trans->get_all_Exons()};
    my $dbID = undef;
    my $rank = 1;
    
    my @exons;
    foreach my $e (@pt_exons) {
      if($e && (!$e->contig || $e->contig->isa('Bio::EnsEMBL::Slice'))) {
	$self->throw('PredictionTranscript must be in contig coords to store');
      }
      
      if($e && $e->isa('Bio::EnsEMBL::StickyExon')) {
	push @exons, @{$e->get_all_component_Exons};
      } else {
	push @exons, $e;
      }
    }

    for my $exon ( @exons ) {
      if( ! defined $exon ) { $rank++; next; }
      
      my $contig_id = $exon->contig->dbID();
      my $contig_start = $exon->start();
      my $contig_end = $exon->end();
      my $contig_strand = $exon->strand();
      my $start_phase = $exon->phase();
      my $end_phase = $exon->end_phase();
      my $score = $exon->score();
      my $p_value = $exon->p_value();
      my $analysis = $pre_trans->analysis->dbID;
      
      if( $rank == 1 ) {
	$exonst->execute( undef, 1, $contig_id, $contig_start, 
			  $contig_end, $contig_strand,
			  $start_phase, $score, $p_value, $analysis, 
			  scalar( @exons ));
	$dbID = $exonst->{'mysql_insertid'};
      } else {
	$exonst->execute( $dbID, $rank, $contig_id, $contig_start, 
			  $contig_end, $contig_strand,
			  $start_phase, $score, $p_value, $analysis, 
			  scalar( @exons ) );
      }
      $rank++;
    }
    
    $pre_trans->dbID( $dbID );
    $pre_trans->adaptor( $self );
  }

  $exonst->finish;
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


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$prediction_transcript_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all prediction transcript features in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("prediction_transcript");
}

1;

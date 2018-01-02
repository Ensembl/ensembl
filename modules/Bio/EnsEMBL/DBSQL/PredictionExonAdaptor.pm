=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::PredictionExonAdaptor - Performs database interaction for
PredictionExons.

=head1 SYNOPSIS

  $pea   = $database_adaptor->get_PredictionExonAdaptor();
  $pexon = $pea->fetch_by_dbID();

  my $slice =
    $database_adaptor->get_SliceAdaptor->fetch_by_region( 'X', 1, 1e6 );

  my @pexons = @{ $pea->fetch_all_by_Slice($slice) };

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::PredictionExonAdaptor;

use vars qw( @ISA );
use strict;


use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::PredictionExon;
use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );


@ISA = qw( Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor );


#_tables
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns the names, aliases of the tables to use for queries
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#

sub _tables {
  return ([ 'prediction_exon', 'pe' ] );
}



#_columns
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns a list of columns to use for queries
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal

sub _columns {
  my $self = shift;

  return qw( pe.prediction_exon_id
             pe.seq_region_id
             pe.seq_region_start
             pe.seq_region_end
             pe.seq_region_strand
             pe.start_phase
             pe.score
             pe.p_value );
}


# _final_clause
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns a default end for the SQL-query (ORDER BY)
#  Returntype : string
#  Exceptions : none
#  Caller     : internal

sub _final_clause {
  return "ORDER BY pe.prediction_transcript_id, pe.exon_rank";
}


=head2 fetch_all_by_PredictionTranscript

  Arg [1]    : Bio::EnsEMBL::PredcitionTranscript $transcript
  Example    : none
  Description: Retrieves all Exons for the Transcript in 5-3 order
  Returntype : listref Bio::EnsEMBL::Exon on Transcript slice 
  Exceptions : throws if transcript does not have a slice
  Caller     : Transcript->get_all_Exons()
  Status     : Stable

=cut

sub fetch_all_by_PredictionTranscript {
  my ( $self, $transcript ) = @_;
  my $constraint = "pe.prediction_transcript_id = ".$transcript->dbID();

  # use 'keep_all' option to keep exons that are off end of slice

  my $tslice = $transcript->slice();
  my $slice;

  if(!$tslice) {
    throw("Transcript must have attached slice to retrieve exons.");
  }

  # use a small slice the same size as the prediction transcript
  $slice = $self->db->get_SliceAdaptor->fetch_by_Feature($transcript);

  my $exons = $self->fetch_all_by_Slice_constraint($slice, $constraint);

  # remap exon coordinates if necessary
  if($slice->name() ne $tslice->name()) {
    my @out;
    foreach my $ex (@$exons) {
      push @out, $ex->transfer($tslice);
    }
    $exons = \@out;
  }

  return $exons;
}



=head2 store

  Arg [1]    : Bio::EnsEMBL::PredictionExon $exon
               The exon to store in this database
  Arg [2]    : int $prediction_transcript_id
               The internal identifier of the prediction exon that that this
               exon is associated with.
  Arg [3]    : int $rank
               The rank of the exon in the transcript (starting at 1)
  Example    : $pexon_adaptor->store($pexon, 1211, 2);
  Description: Stores a PredictionExon in the database
  Returntype : none
  Exceptions : thrown if exon does not have a slice attached
               or if $exon->start, $exon->end, $exon->strand, or $exon->phase 
               are not defined or if $exon is not a Bio::EnsEMBL::PredictionExon 
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ( $self, $pexon, $pt_id, $rank ) = @_;

  if(!ref($pexon) || !$pexon->isa('Bio::EnsEMBL::PredictionExon') ) {
    throw("Expected PredictionExon argument");
  }

  throw("Expected PredictionTranscript id argument.") if(!$pt_id);
  throw("Expected rank argument.") if(!$rank);

  my $db = $self->db();

  if($pexon->is_stored($db)) {
    warning('PredictionExon is already stored in this DB.');
    return $pexon->dbID();
  }

  if( ! $pexon->start || ! $pexon->end ||
      ! $pexon->strand || ! defined $pexon->phase ) {
    throw("PredictionExon does not have all attributes to store.\n" .
         "start, end, strand and phase attributes must be set.");
  }

  #maintain reference to original passed-in prediction exon
  my $original = $pexon;
  my $seq_region_id;
  ($pexon, $seq_region_id) = $self->_pre_store($pexon);

  my $sth = $db->dbc->prepare
    ("INSERT into prediction_exon (prediction_transcript_id, exon_rank, " .
                       "seq_region_id, seq_region_start, seq_region_end, " .
                       "seq_region_strand, start_phase, score, p_value) " .
      "VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ? )");

  $sth->bind_param(1,$pt_id,SQL_INTEGER);
  $sth->bind_param(2,$rank,SQL_SMALLINT);
  $sth->bind_param(3,$seq_region_id,SQL_INTEGER);
  $sth->bind_param(4,$pexon->start,SQL_INTEGER);
  $sth->bind_param(5,$pexon->end,SQL_INTEGER);
  $sth->bind_param(6,$pexon->strand,SQL_TINYINT);
  $sth->bind_param(7,$pexon->phase,SQL_TINYINT);
  $sth->bind_param(8,$pexon->score,SQL_DOUBLE);
  $sth->bind_param(9,$pexon->p_value,SQL_DOUBLE);

  $sth->execute();

  my $dbID = $self->last_insert_id('prediction_transcript_id', undef, 'prediction_transcript');

  #set the adaptor and dbID of the object they passed in
  $original->dbID($dbID);
  $original->adaptor($self);

  return $dbID;
}



=head2 remove

  Arg [1]    : Bio::EnsEMBL::PredictionExon $exon
               the exon to remove from the database 
  Example    : $exon_adaptor->remove($exon);
  Description: Removes an exon from the database
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $pexon = shift;

  my $db = $self->db();

  if(!$pexon->is_stored($db)) {
    warning('PredictionExon is not in this DB - not removing');
    return undef;
  }

  my $sth = $self->prepare(
            "DELETE FROM prediction_exon WHERE prediction_exon_id = ?");
  $sth->bind_param( 1, $pexon->dbID, SQL_INTEGER );
  $sth->execute();

  $pexon->dbID(undef);
  $pexon->adaptor(undef);
}



=head2 list_dbIDs

  Arg [1]    : none
  Example    : @exon_ids = @{$exon_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all exons in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self,$ordered) = @_;

   return $self->_list_dbIDs("prediction_exon",undef, $ordered);
}



#_objs_from_sth

#  Arg [1]    : Hashreference $hashref
#  Example    : none 
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of Genes 
#  Returntype : listref of Bio::EnsEMBL::Genes in target coordinate system
#  Exceptions : none
#  Caller     : internal
#

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #
  my $sa = $self->db()->get_SliceAdaptor();

  my @exons;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my(
    $prediction_exon_id, $seq_region_id,
    $seq_region_start,   $seq_region_end, $seq_region_strand,
    $start_phase,        $score,          $p_value);

  $sth->bind_columns(\(
                      $prediction_exon_id, $seq_region_id,
                      $seq_region_start,   $seq_region_end, $seq_region_strand,
                      $start_phase,        $score,          $p_value));

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  my $dest_slice_cs;
  my $dest_slice_sr_name;
  my $dest_slice_sr_id;
  my $asma;
 
  if ($dest_slice) {
    $dest_slice_start   = $dest_slice->start();
    $dest_slice_end     = $dest_slice->end();
    $dest_slice_strand  = $dest_slice->strand();
    $dest_slice_length  = $dest_slice->length();
    $dest_slice_cs      = $dest_slice->coord_system();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
    $dest_slice_sr_id   = $dest_slice->get_seq_region_id();
    $asma               = $self->db->get_AssemblyMapperAdaptor();
  }

  FEATURE: while($sth->fetch()) {

    #need to get the internal_seq_region, if present
    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
    my $slice = $slice_hash{"ID:".$seq_region_id};

    if (!$slice) {
      $slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id}       = $slice->coord_system();
    }

    #obtain a mapper if none was defined, but a dest_seq_region was
    if(!$mapper && $dest_slice && !$dest_slice_cs->equals($slice->coord_system)) {
      $mapper = $asma->fetch_by_CoordSystems($dest_slice_cs, $slice->coord_system);
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};

    #
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    #

    if ($mapper) {

      if (defined $dest_slice && $mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper') ) {
        ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) =
         $mapper->map($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs, 1, $dest_slice);

      } else {
        ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) =
         $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);
      }

      #skip features that map to gaps or coord system boundaries
      next FEATURE if (!defined($seq_region_id));

      #get a slice in the coord system we just mapped to
      $slice = $slice_hash{"ID:".$seq_region_id} ||= $sa->fetch_by_seq_region_id($seq_region_id);
    }

    #
    # If a destination slice was provided convert the coords.
    #
    if (defined($dest_slice)) {
      my $seq_region_len = $dest_slice->seq_region_length();

      if ( $dest_slice_strand == 1 ) {
        $seq_region_start = $seq_region_start - $dest_slice_start + 1;
        $seq_region_end   = $seq_region_end - $dest_slice_start + 1;

        if ( $dest_slice->is_circular ) {
        # Handle circular chromosomes.

          if ( $seq_region_start > $seq_region_end ) {
            # Looking at a feature overlapping the chromosome origin.

            if ( $seq_region_end > $dest_slice_start ) {
              # Looking at the region in the beginning of the chromosome
              $seq_region_start -= $seq_region_len;
            }
            if ( $seq_region_end < 0 ) {
              $seq_region_end += $seq_region_len;
            }
          } else {
            if ($dest_slice_start > $dest_slice_end && $seq_region_end < 0) {
              # Looking at the region overlapping the chromosome
              # origin and a feature which is at the beginning of the
              # chromosome.
              $seq_region_start += $seq_region_len;
              $seq_region_end   += $seq_region_len;
            }
          }
        }
      } else {

        my $start = $dest_slice_end - $seq_region_end + 1;
        my $end = $dest_slice_end - $seq_region_start + 1;

        if ($dest_slice->is_circular()) {

          if ($dest_slice_start > $dest_slice_end) {
            # slice spans origin or replication

            if ($seq_region_start >= $dest_slice_start) {
              $end += $seq_region_len;
              $start += $seq_region_len if $seq_region_end > $dest_slice_start;

            } elsif ($seq_region_start <= $dest_slice_end) {
              # do nothing
            } elsif ($seq_region_end >= $dest_slice_start) {
              $start += $seq_region_len;
              $end += $seq_region_len;

            } elsif ($seq_region_end <= $dest_slice_end) {
              $end += $seq_region_len if $end < 0;

            } elsif ($seq_region_start > $seq_region_end) {
              $end += $seq_region_len;
            }

          } else {

            if ($seq_region_start <= $dest_slice_end and $seq_region_end >= $dest_slice_start) {
              # do nothing
            } elsif ($seq_region_start > $seq_region_end) {
              if ($seq_region_start <= $dest_slice_end) {
                $start -= $seq_region_len;
              } elsif ($seq_region_end >= $dest_slice_start) {
                $end += $seq_region_len;
              }
            }
          }
        }

        $seq_region_start = $start;
        $seq_region_end = $end;
        $seq_region_strand *= -1;

      } ## end else [ if ( $dest_slice_strand...)]

      # Throw away features off the end of the requested slice or on
      # different seq_region.
      if ($seq_region_end < 1
          || $seq_region_start > $dest_slice_length
          || ($dest_slice_sr_id != $seq_region_id)) {
        next FEATURE;
      }
      $slice = $dest_slice;
    }

    # Finally, create the new PredictionExon.
    push( @exons,
          $self->_create_feature( 'Bio::EnsEMBL::PredictionExon', {
                                    '-start'   => $seq_region_start,
                                    '-end'     => $seq_region_end,
                                    '-strand'  => $seq_region_strand,
                                    '-adaptor' => $self,
                                    '-slice'   => $slice,
                                    '-dbID'    => $prediction_exon_id,
                                    '-phase'   => $start_phase,
                                    '-score'   => $score,
                                    '-p_value' => $p_value
                                  } ) );

  }

  return \@exons;
}


1;

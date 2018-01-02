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

Bio::EnsEMBL::Map::DBSQL::DitagFeatureAdaptor

=head1 SYNOPSIS

  my $dfa = $db->get_DitagFeatureAdaptor;
  my $ditagFeatures = $dfa->fetch_all_by_Slice( $slice, "SME005" );

  foreach my $ditagFeature (@$ditagFeatures) {
    print $ditagFeature->ditag_id . " "
      . $ditagFeature->slice . " "
      . $ditagFeature->start . "-"
      . $ditagFeature->end . " "
      . $ditagFeature->strand;
  }

=head1 DESCRIPTION

Provides database interaction for the Bio::EnsEMBL::Map::DitagFeature
object

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::DBSQL::DitagFeatureAdaptor;

use strict;
use vars ('@ISA');

use Bio::EnsEMBL::Map::Ditag;
use Bio::EnsEMBL::Map::DitagFeature;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_all

  Arg [1]    : none
  Example    : @all_tags = @{$ditagfeature_adaptor->fetch_all};
  Description: Retrieves all ditagFeatures from the database;
               Usually not a good idea, use fetch_all_by_Slice instead.
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeature
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all {
  my $self = shift;

  my $sth = $self->prepare("SELECT df.ditag_feature_id, df.ditag_id, df.seq_region_id, 
                            df.seq_region_start, df.seq_region_end, df.seq_region_strand, 
                            df.analysis_id, df.hit_start, df.hit_end, df.hit_strand, 
                            df.cigar_line, df.ditag_side, df.ditag_pair_id, d.tag_count 
                            FROM   ditag_feature df, ditag d 
                            WHERE  df.ditag_id=d.ditag_id" );
  $sth->execute;

  my $result = $self->_fetch($sth);

  return $result;
}


=head2 fetch_by_dbID

  Arg [1]    : ditagFeature dbID
  Example    : @my_tags = @{$ditagfeature_adaptor->fetch_by_dbID($my_id)};
  Description: Retrieves a ditagFeature from the database.
  Returntype : Bio::EnsEMBL::Map::DitagFeature
  Caller     : general
  Status     : At Risk

=cut

sub fetch_by_dbID {
  my ($self, $dbid) = @_;

  my $sth = $self->prepare("SELECT df.ditag_feature_id, df.ditag_id, df.seq_region_id, 
                            df.seq_region_start, df.seq_region_end, df.seq_region_strand, 
                            df.analysis_id, df.hit_start, df.hit_end, df.hit_strand, 
                            df.cigar_line, df.ditag_side, df.ditag_pair_id, d.tag_count 
                            FROM   ditag_feature df, ditag d 
                            WHERE  df.ditag_id=d.ditag_id AND df.ditag_feature_id = ?" );
  $sth->execute($dbid);

  my $result = $self->_fetch($sth);

  return $result->[0];
}


=head2 fetch_all_by_ditagID

  Arg [1]    : ditag dbID
  Arg [2]    : (optional) ditag-pair dbID
  Arg [3]    : (optional) analysis ID
  Example    : @my_tags = @{$ditagfeature_adaptor->fetch_all_by_ditag_id($my_id)};
  Description: Retrieves all ditagFeatures from the database linking to a specific ditag-id
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeature
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_ditagID {
  my ($self, $ditag_id, $ditag_pair_id, $analysis_id) = @_;

  my $arg = $ditag_id;
  my $sql = "SELECT df.ditag_feature_id, df.ditag_id, df.seq_region_id, 
             df.seq_region_start, df.seq_region_end, df.seq_region_strand, 
             df.analysis_id, df.hit_start, df.hit_end, df.hit_strand, 
             df.cigar_line, df.ditag_side, df.ditag_pair_id, d.tag_count 
             FROM   ditag_feature df, ditag d 
             WHERE  df.ditag_id=d.ditag_id AND df.ditag_id = ? ";
  if($ditag_pair_id){
    $sql .= "AND df.ditag_pair_id = ? ";
    $arg .= ", $ditag_pair_id";
  }
  if($analysis_id){
    $sql .= "AND df.analysis_id = ? ";
    $arg .= ", $analysis_id";
  }
  $sql   .= "ORDER BY df.ditag_pair_id";
  my $sth = $self->prepare($sql);
  $sth->execute(split(",",$arg));

  my $result = $self->_fetch($sth);

  return $result;
}


=head2 fetch_all_by_type

  Arg [1]    : ditag type
  Example    : @my_tags = @{$ditagfeature_adaptor->fetch_all_by_type($type)};
  Description: Retrieves all ditagFeatures from the database linking to a specific ditag-type
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeature
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_type {
  my ($self, $ditag_type) = @_;

  my $sth = $self->prepare("SELECT df.ditag_feature_id, df.ditag_id, df.seq_region_id, 
                            df.seq_region_start, df.seq_region_end, df.seq_region_strand, 
                            df.analysis_id, df.hit_start, df.hit_end, df.hit_strand, 
                            df.cigar_line, df.ditag_side, df.ditag_pair_id, d.tag_count 
                            FROM   ditag_feature df, ditag d 
                            WHERE  df.ditag_id=d.ditag_id AND d.type = ? 
                            ORDER BY df.ditag_id, df.ditag_pair_id" );
  $sth->execute($ditag_type);

  my $result = $self->_fetch($sth);

  return $result;
}


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : (optional) ditag type name (specific library) or an aray ref with multiple type names
  Arg [3]    : (optional) analysis logic_name
  Example    : $tags = $ditagfeature_adaptor->fetch_all_by_Slice($slice, "SME005");
  Description: Retrieves ditagFeatures from the database overlapping a specific region
               and (optional) of a specific ditag type or analysis.
               Start & end locations are returned in slice coordinates, now.
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeatures
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $tagtype, $logic_name) = @_;

  my @result;
  my $moresql;

  if(!ref($slice) || !$slice->isa("Bio::EnsEMBL::Slice")) {
    throw("Bio::EnsEMBL::Slice argument expected not $slice.");
  }

  #get affected ditag_feature_ids
  my $sql = "SELECT df.ditag_feature_id, df.ditag_id, df.seq_region_id, df.seq_region_start, 
             df.seq_region_end, df.seq_region_strand, df.analysis_id, df.hit_start, df.hit_end, 
             df.hit_strand, df.cigar_line, df.ditag_side, df.ditag_pair_id, 
             d.tag_count 
             FROM ditag_feature df, ditag d 
             WHERE df.ditag_id=d.ditag_id";
  if($tagtype){
    my $tagtypes = '';
    #check if array
    if(ref $tagtype eq 'ARRAY'){
      my @arraytype_mod;
      foreach my $arraytype (@$tagtype){ push @arraytype_mod, '"'.$arraytype.'"' }
      $tagtypes = join(", ", @arraytype_mod);
    }
    else{
      $tagtypes = '"'.$tagtype.'"';
    }
    $sql .= " AND d.type IN(".$tagtypes.")";
  }
  if($logic_name){
    my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
    if(!$analysis) {
      return undef;
    }
    $sql .= " AND df.analysis_id = ".$analysis->dbID();
  }
  $sql .= " AND df.seq_region_id = ".$slice->get_seq_region_id.
          " AND df.seq_region_start <= ".$slice->end.
	  " AND df.seq_region_end   >= ".$slice->start;

  my $sth = $self->prepare($sql);
  $sth->execute();

  my $result = $self->_fetch($sth, $slice);
  push(@result, @$result);

  return \@result;
}


=head2 fetch_pairs_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : (optional) ditag type (specific library)
  Arg [3]    : (optional) analysis logic_name
  Example    : my $ditagfeatures = $dfa->fetch_pairs_by_Slice($slice);
               foreach my $ditagfeature (@$ditagfeatures){
                 $minstart   = $$ditagfeature2{'start'};
                 $maxend     = $$ditagfeature2{'end'};
                 $bothstrand = $$ditagfeature2{'strand'};
                 $tag_count  = $$ditagfeature2{'tag_count'};
                 print "$minstart, $maxend, $bothstrand, $tag_count\n";
               }
  Description: Retrieves ditagFeature information in pairs from the database overlapping a specific region
               and (optional) of a specific ditag type or analysis. The absotute start and end points are
               fetched.
               Slices should be SMALL!
  Returntype : array ref with hash ref of artifical DitagFeature object
  Caller     : general
  Status     : At Risk

=cut

sub fetch_pairs_by_Slice {
  my ($self, $slice, $tagtype, $logic_name) = @_;
  my ($tag_id, $pair_id, $seq_region_id, $start, $end, $strand, $analysis_id, $tag_count);
  my @result;

  my $sql = "SELECT df.ditag_id, df.ditag_pair_id, df.seq_region_id, MIN(df.seq_region_start), ".
            "MAX(df.seq_region_end), df.seq_region_strand, df.analysis_id, d.tag_count ".
            "FROM ditag_feature df, ditag d ".
            "WHERE df.ditag_id=d.ditag_id ";
  if($tagtype){
    $sql .= "AND d.type = \"".$tagtype."\"";
  }
  $sql .= " AND df.seq_region_id = ".$slice->get_seq_region_id.
          " AND df.seq_region_start <= ".$slice->end.
	  " AND df.seq_region_end >= ".$slice->start;
  if($logic_name){
    my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
    if(!$analysis) {
      return undef;
    }
    $sql .= " AND df.analysis_id = ".$analysis->dbID();
  }
  $sql .= " GROUP BY df.ditag_id, df.ditag_pair_id;";
  my $sth = $self->prepare($sql);
  $sth->execute();
  $sth->bind_columns( \$tag_id, \$pair_id, \$seq_region_id, \$start, \$end, \$strand, \$analysis_id ,\$tag_count);
  while ( $sth->fetch ) {
    # convert into relative slice coordinates
    my $seq_region_len = $slice->seq_region_length();

    if ($slice->strand == 1) { # Positive strand
      $start = $start - $slice->start + 1;
      $end   = $end - $slice->start + 1;

      if ($slice->is_circular()) { # Handle circular chromosomes
	if ($start > $end) { # Looking at a feature overlapping the chromsome origin
	  if ($end > $slice->start) {
	    # Looking at the region in the beginning of the chromosome
	    $start -= $seq_region_len;
	  }

	  if ($end < 0) {
	    $end += $seq_region_len;
	  }
	} else {
	  if ($slice->start > $slice->end && $end < 0) {
	    # Looking at the region overlapping the chromosome origin and 
	    # a feature which is at the beginning of the chromosome
	    $start += $seq_region_len;
	    $end   += $seq_region_len;
	  }
	}
      } # end if ($dest_slice->is_circular...)

    } else { # Negative strand
      my ($seq_region_start, $seq_region_end) = ($start, $end);
      $start = $slice->end - $seq_region_end + 1;
      $end = $slice->end - $seq_region_start + 1;

      if ($slice->is_circular()) {
	if ($slice->start > $slice->end) { # slice spans origin or replication
	  if ($seq_region_start >= $slice->start) {
	    $end += $seq_region_len;
	    $start += $seq_region_len 
	      if $seq_region_end > $slice->start;

	  } elsif ($seq_region_start <= $slice->end) {
	    # do nothing
	  } elsif ($seq_region_end >= $slice->start) {
	    $start += $seq_region_len;
	    $end += $seq_region_len;

	  } elsif ($seq_region_end <= $slice->end) {
	    $end += $seq_region_len
	      if $end < 0;
	  } elsif ($seq_region_start > $seq_region_end) {
	    $end += $seq_region_len;
	  } else { }
      
	} else {
	  if ($seq_region_start <= $slice->end and $seq_region_end >= $slice->start) {
	    # do nothing
	  } elsif ($seq_region_start > $seq_region_end) {
	    if ($seq_region_start <= $slice->end) {
	      $start -= $seq_region_len;
	    } elsif ($seq_region_end >= $slice->start) {
	      $end += $seq_region_len;
	    } else { }
	  }
	}
      }

      $strand *= -1;
    }

    my %ditag_feature_pair = (
                      ditag     => $tag_id,
                      pair_id   => $pair_id,
                      region    => $seq_region_id,
                      start     => $start,
                      end       => $end,
                      strand    => $strand,
                      analysis  => $analysis_id,
                      tag_count => $tag_count
                     );
    push(@result, \%ditag_feature_pair);
  }

  return \@result;
}


=head2 _fetch

  Arg [1]    : statement handler
  Arg [2]    : (optional) target-slice for the feature
  Description: generic sql-fetch function for the DitagFeature fetch methods
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeatures
  Caller     : private
  Status     : At Risk

=cut

sub _fetch {
  my ($self, $sth, $dest_slice) = @_;

  my ( $tag_id, $mothertag_id, $seqreg, $seqstart, $seqend, $strand, $analysis_id, $hit_start,
       $hit_end, $hit_strand, $cigar_line, $ditag_side, $ditag_pair_id, $tag_count );
  $sth->bind_columns( \$tag_id,        \$mothertag_id, \$seqreg,
                      \$seqstart,      \$seqend,       \$strand,
                      \$analysis_id,   \$hit_start,    \$hit_end,
                      \$hit_strand,    \$cigar_line,   \$ditag_side,
                      \$ditag_pair_id, \$tag_count );

  my @ditag_features;
  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  if($dest_slice) {
    $dest_slice_start   = $dest_slice->start();
    $dest_slice_end     = $dest_slice->end();
    $dest_slice_strand  = $dest_slice->strand();
  }

  while ( $sth->fetch ) {
    my $analysis_obj = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
    my $slice        = $self->db->get_SliceAdaptor->fetch_by_seq_region_id($seqreg);

    if($dest_slice) {
      if($dest_slice_start != 1 || $dest_slice_strand != 1) {
        if($dest_slice_strand == 1) {
          $seqstart    = $seqstart  - $dest_slice_start + 1;
          $seqend      = $seqend    - $dest_slice_start + 1;
        } else {
          my $tmp_seq_region_start = $seqstart;
          $seqstart    = $dest_slice_end - $seqend + 1;
          $seqend      = $dest_slice_end - $tmp_seq_region_start + 1;
          $strand     *= -1;
        }
      }

      my $seq_region_len = $dest_slice->seq_region_length();

      if ($dest_slice_strand == 1) { # Positive strand		
	$seqstart = $seqstart - $dest_slice_start + 1;
	$seqend   = $seqend - $dest_slice_start + 1;

	if ($dest_slice->is_circular()) { # Handle cicular chromosomes
	  if ($seqstart > $seqend) { # Looking at a feature overlapping the chromsome origin
	    if ($seqend > $dest_slice_start) {
	      # Looking at the region in the beginning of the chromosome.
	      $seqstart -= $seq_region_len;
	    }

	    if ($seqend < 0) {
	      $seqend += $seq_region_len;
	    }
	  } else {
	    if ($dest_slice_start > $dest_slice_end && $seqend < 0) {
	      # Looking at the region overlapping the chromosome origin and 
	      # a feature which is at the beginning of the chromosome.
	      $seqstart += $seq_region_len;
	      $seqend   += $seq_region_len;
	    }
	  }
	}
      } else { # Negative strand
	my $start = $dest_slice_end - $seqend + 1;
	my $end = $dest_slice_end - $seqstart + 1;

	if ($dest_slice->is_circular()) {
	  if ($dest_slice_start > $dest_slice_end) { 
	    # slice spans origin or replication
	    if ($seqstart >= $dest_slice_start) {
	      $end += $seq_region_len;
	      $start += $seq_region_len 
		if $seqend > $dest_slice_start;

	    } elsif ($seqstart <= $dest_slice_end) {
	      # do nothing
	    } elsif ($seqend >= $dest_slice_start) {
	      $start += $seq_region_len;
	      $end += $seq_region_len;
	    } elsif ($seqend <= $dest_slice_end) {
	      $end += $seq_region_len
		if $end < 0;
	    } elsif ($seqstart > $seqend) {
	      $end += $seq_region_len;
	    } else { }
	  } else {
	    if ($seqstart <= $dest_slice_end and $seqend >= $dest_slice_start) {
	      # do nothing
	    } elsif ($seqstart > $seqend) {
	      if ($seqstart <= $dest_slice_end) {
		$start -= $seq_region_len;
	      } elsif ($seqend >= $dest_slice_start) {
		$end += $seq_region_len;
	      } else { }
	    }
	  }
	}

	$seqstart = $start;
	$seqend = $end;
	$strand *= -1;
      }

      $slice = $dest_slice;
    }

    push @ditag_features,
      Bio::EnsEMBL::Map::DitagFeature->new( -dbid          => $tag_id,
                                            -slice         => $slice,
                                            -start         => $seqstart,
                                            -end           => $seqend,
                                            -strand        => $strand, 
                                            -analysis      => $analysis_obj,
                                            -hit_start     => $hit_start,
                                            -hit_end       => $hit_end,
                                            -hit_strand    => $hit_strand,
                                            -ditag_id      => $mothertag_id,
                                            -cigar_line    => $cigar_line,
                                            -ditag_side    => $ditag_side,
					    -ditag_pair_id => $ditag_pair_id,
					    -ditag         => undef,
					    -tag_count     => $tag_count,
                                            -adaptor       => $self,
                                            );
  }

  return \@ditag_features;
}


=head2 sequence

  Arg [1]    : dbID of DitagFeature
  Example    : $ditagfeature_adaptor->get_sequence($ditagFeature->dbID)
  Description: get the part of the sequence of a ditag,
               that is actully aligned to the genome.
  Returntype : string
  Exceptions : thrown if not all data needed for storing is populated in the
               ditag features
  Caller     : Bio::EnsEMBL::Map::DitagFeature
  Status     : At Risk

=cut

sub sequence {
  my ($self, $dbID) = @_;

  my $sequence = undef;
  my $db  = $self->db() or throw "Couldn t get database connection.";
  my $sql = "SELECT d.sequence, df.hit_start, df.hit_end, df.hit_strand ".
            "FROM ditag d, ditag_feature df ".
	    "WHERE df.ditag_id=d.ditag_id and df.ditag_feature_id = ?";
  my $sth = $db->dbc->prepare($sql);
  $sth->execute( $dbID );
  my ($seq, $start, $end, $strand) = $sth->fetchrow_array();
  if($seq and $start and $end and $strand){
    $sequence = substr($seq, ($start-1), ($end-$strand));
    if($strand == -1) {
      $sequence =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    }
  }

  return $sequence;
}


=head2 store

  Arg [1]    : (Array ref of) Bio::EnsEMBL::Map::DitagFeature
  Example    : $ditagfeature_adaptor->store(@ditag_features);
  Description: Stores a single ditagFeature or
               a list of ditagFeatures in this database.
  Returntype : none
  Exceptions : thrown if not all data needed for storing is populated in the
               ditag features
  Caller     : general
  Status     : At Risk

=cut

sub store {
  my ( $self, $ditag_features ) = @_;

  if ( ref $ditag_features eq 'ARRAY' ) {
    if ( scalar(@$ditag_features) == 0 ) {
      throw( "Must call store with ditag_feature or list ref of ditags_features" );
    }
  } elsif ($ditag_features) {
    my @ditag_features;
    push @ditag_features, $ditag_features;
    $ditag_features = \@ditag_features;
  } else {
    throw( "Must call store with ditag_feature or list ref of ditag_features." );
  }

  my $db = $self->db() or throw "Couldn t get database connection.";

  my $sth1 = $self->prepare( "INSERT INTO ditag_feature( ditag_id, seq_region_id, seq_region_start, 
                              seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, 
                              hit_strand, cigar_line, ditag_side, ditag_pair_id ) 
                              VALUES( ?,?,?,?,?,?,?,?,?,?,?,? )" );
  my $sth2 = $self->prepare( "INSERT INTO ditag_feature( ditag_feature_ID, ditag_id, seq_region_id, 
                              seq_region_start, seq_region_end, seq_region_strand, analysis_id, hit_start, 
                              hit_end, hit_strand, cigar_line, ditag_side, ditag_pair_id ) 
                              VALUES( ?,?,?,?,?,?,?,?,?,?,?,?,? )" );
#  my $sth3 = $self->prepare( "SELECT COUNT(*) FROM ditag_feature 
#                              WHERE ditag_id = ?" );

TAG:
  foreach my $ditag_feature (@$ditag_features) {

    if ( !ref $ditag_feature || !$ditag_feature->isa("Bio::EnsEMBL::Map::DitagFeature") ) {
      throw(   "Object must be an Ensembl DitagFeature, "
             . "not a " . ref($ditag_feature) );
    }
    if ( $ditag_feature->is_stored($db) ) {
      warning(   "DitagFeature " . $ditag_feature->dbID .
                 " is already stored in this database,".
                 " maybe you need to use the update() method?" );
      next TAG;
    }
    if(!$ditag_feature->ditag_id or !($self->db->get_DitagAdaptor->fetch_by_dbID($ditag_feature->ditag_id))){
      throw("DitagFeature must be supplied with the id of a corresponding Ditag object.");
    }
    if(!$ditag_feature->ditag or !$ditag_feature->ditag->isa("Bio::EnsEMBL::Map::Ditag")){
      throw("DitagFeature must be linked to a valid Ditag object.");
    }


#    #check if more than x tags with this ditag id exist
#    $sth3->execute( $ditag_feature->ditag_id );
#    my ($num) = $sth3->fetchrow_array();
#    if ( ($num) and ($num > 1) ) {
#      warning( "There are already at least 2 DitagFeatures relating to Ditag ".
#                 $ditag->ditag_id." stored in this database." );
#      if ( $num > 4 ) {
#        warning( "not storing" );
#        next TAG;
#      }
#    }

    if ( $ditag_feature->dbID ) {
      $sth2->bind_param( 1,  $ditag_feature->dbID,                      SQL_INTEGER );
      $sth2->bind_param( 2,  $ditag_feature->ditag_id,                  SQL_INTEGER );
      $sth2->bind_param( 3, ($ditag_feature->slice->get_seq_region_id), SQL_INTEGER );
      $sth2->bind_param( 4,  $ditag_feature->start,                     SQL_INTEGER );
      $sth2->bind_param( 5,  $ditag_feature->end,                       SQL_INTEGER );
      $sth2->bind_param( 6,  $ditag_feature->strand,                    SQL_VARCHAR );
      $sth2->bind_param( 7,  $ditag_feature->analysis->dbID,            SQL_INTEGER );
      $sth2->bind_param( 8,  $ditag_feature->hit_start,                 SQL_INTEGER );
      $sth2->bind_param( 9,  $ditag_feature->hit_end,                   SQL_INTEGER );
      $sth2->bind_param( 10, $ditag_feature->hit_strand,                SQL_VARCHAR );
      $sth2->bind_param( 11, $ditag_feature->cigar_line,                SQL_VARCHAR );
      $sth2->bind_param( 12, $ditag_feature->ditag_side,                SQL_VARCHAR );
      $sth2->bind_param( 13, $ditag_feature->ditag_pair_id,             SQL_VARCHAR );
      $sth2->execute();
    }
    else{
      $sth1->bind_param( 1,  $ditag_feature->ditag_id,                  SQL_INTEGER );
      $sth1->bind_param( 2, ($ditag_feature->slice->get_seq_region_id), SQL_INTEGER );
      $sth1->bind_param( 3,  $ditag_feature->start,                     SQL_INTEGER );
      $sth1->bind_param( 4,  $ditag_feature->end,                       SQL_INTEGER );
      $sth1->bind_param( 5,  $ditag_feature->strand,                    SQL_VARCHAR );
      $sth1->bind_param( 6,  $ditag_feature->analysis->dbID,            SQL_INTEGER );
      $sth1->bind_param( 7,  $ditag_feature->hit_start,                 SQL_INTEGER );
      $sth1->bind_param( 8,  $ditag_feature->hit_end,                   SQL_INTEGER );
      $sth1->bind_param( 9,  $ditag_feature->hit_strand,                SQL_VARCHAR );
      $sth1->bind_param( 10, $ditag_feature->cigar_line,                SQL_VARCHAR );
      $sth1->bind_param( 11, $ditag_feature->ditag_side,                SQL_VARCHAR );
      $sth1->bind_param( 12, $ditag_feature->ditag_pair_id,             SQL_VARCHAR );
      $sth1->execute();
      my $dbID = $self->last_insert_id('ditag_feature_id', undef, 'ditag_feature');
      $ditag_feature->dbID($dbID);
      $ditag_feature->adaptor($self);
    }

  }
}


=head2 batch_store

  Arg [1]    : (Array ref of) Bio::EnsEMBL::Map::DitagFeatures
  Arg [2]    : bool have_dbIDs
  Example    : $ditagfeature_adaptor->batch_store(\@ditag_features);
  Description: Stores a list of ditagFeatures in this database.
               DitagFeatures are expected to have no dbID yet unless flag "have_dbIDs" is true.
               They are inserted in one combined INSERT for better performance.
  Returntype : none
  Exceptions : thrown if not all data needed for storing is given for the
               ditag features
  Caller     : general
  Status     : At Risk

=cut

sub batch_store {
  my ( $self, $ditag_features, $have_dbIDs ) = @_;

  my @good_ditag_features;
  my ($sql, $sqladd);
  my $inserts = 0;

  if ( ref $ditag_features eq 'ARRAY' ) {
    if ( scalar(@$ditag_features) == 0 ) {
      throw( "Must call store with ditag_feature or list ref of ditag_features." );
    }
  } elsif ($ditag_features) {
    my @ditag_features;
    push @ditag_features, $ditag_features;
    $ditag_features = \@ditag_features;
  } else {
    throw( "Must call store with ditag_feature or list ref of ditag_features." );
  }

  my $db = $self->db() or throw "Couldn t get database connection.";

  #check whether it s a DitagFeature object and is not stored already
  foreach my $ditag_feature (@$ditag_features) {

    if ( !ref $ditag_feature || !$ditag_feature->isa("Bio::EnsEMBL::Map::DitagFeature") ) {
      throw(   "Object must be an Ensembl DitagFeature, "
             . "not a " . ref($ditag_feature) );
    }
    if(!$ditag_feature->ditag_id or !($self->db->get_DitagAdaptor->fetch_by_dbID($ditag_feature->ditag_id))){
      throw("DitagFeature must be supplied with the id of a corresponding Ditag object.");
    }

    if(!$ditag_feature->ditag or !$ditag_feature->ditag->isa("Bio::EnsEMBL::Map::Ditag")){
      throw("DitagFeature must be linked to a valid Ditag object.");
    }
    if ( $ditag_feature->is_stored($db) ) {
      warning(   "DitagFeature " . $ditag_feature->dbID
                 . " is already stored in this database." );
      next;
    }
    push(@good_ditag_features, $ditag_feature);
  }
  $ditag_features = undef;

  #create batch INSERT
  if($have_dbIDs){
    $sql = "INSERT INTO ditag_feature ( ditag_feature_id, ditag_id, seq_region_id, seq_region_start, ".
           "seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, ".
           "hit_strand, cigar_line, ditag_side, ditag_pair_id ) VALUES ";
    foreach my $ditag_feature (@good_ditag_features) {
      $sqladd = "";
      if($inserts){ $sqladd = ", " }
      $sqladd .= "(". $ditag_feature->ditag_feature_id.", ".$ditag_feature->ditag_id.", ".
	         ($ditag_feature->slice->get_seq_region_id).", ". $ditag_feature->start.", ".
		 $ditag_feature->end.", '".$ditag_feature->strand."', ".$ditag_feature->analysis->dbID.", ".
                 $ditag_feature->hit_start.", ".$ditag_feature->hit_end.", '".$ditag_feature->hit_strand.
		 "', '".$ditag_feature->cigar_line."', '".$ditag_feature->ditag_side."', ".
		 $ditag_feature->ditag_pair_id.")";
      $sql .= $sqladd;
      $inserts++;
    }
  }
  else{
    $sql = "INSERT INTO ditag_feature ( ditag_id, seq_region_id, seq_region_start, ".
           "seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, ".
           "hit_strand, cigar_line, ditag_side, ditag_pair_id ) VALUES ";
    foreach my $ditag_feature (@good_ditag_features) {
      $sqladd = "";
      if($inserts){ $sqladd = ", " }
      $sqladd .= "(". $ditag_feature->ditag_id.", ".($ditag_feature->slice->get_seq_region_id).", ".
	         $ditag_feature->start.", ".$ditag_feature->end.", '".$ditag_feature->strand."', ".
		 $ditag_feature->analysis->dbID.", ".$ditag_feature->hit_start.", ".$ditag_feature->hit_end.
		 ", '".$ditag_feature->hit_strand."', '".$ditag_feature->cigar_line."', '".
		 $ditag_feature->ditag_side."', ".$ditag_feature->ditag_pair_id.")";
      $sql .= $sqladd;
      $inserts++;
    }
  }

  #STORE
  if($inserts){
    print STDERR "\nHave $inserts Features.\n";
    eval{
      $db->dbc->do($sql);
    };
    if($@){
      warning("Problem inserting ditag feature batch!".$@."\n");
    }
  }
  else{
    warn "Nothing stored!";
  }

}


=head2 update

  Arg [1]    : ditagFeature to update
  Description: update an existing ditagFeature with new values
  Returntype : 1 on success
  Status     : At Risk

=cut

sub update {
  my ($self, $ditagFeature) = @_;

  my $sth = $self->prepare( "UPDATE ditag_feature
                             SET ditag_id=?, seq_region_id=?, seq_region_start=?, seq_region_end=?,
                             seq_region_strand=?, analysis_id=?, hit_start=?, hit_end=?, hit_strand=?,
                             cigar_line=?, ditag_side=?, ditag_pair_id=?
                             where ditag_feature_id=?;" );

  $sth->bind_param(1, $ditagFeature->ditag_id, SQL_INTEGER);
  $sth->bind_param(1, $ditagFeature->seq_region_id, SQL_INTEGER);
  $sth->bind_param(1, $ditagFeature->seq_region_start, SQL_INTEGER);
  $sth->bind_param(1, $ditagFeature->seq_region_end, SQL_INTEGER);
  $sth->bind_param(1, $ditagFeature->seq_region_strand, SQL_TINYINT);
  $sth->bind_param(1, $ditagFeature->hit_start, SQL_INTEGER);
  $sth->bind_param(1, $ditagFeature->hit_end, SQL_INTEGER);
  $sth->bind_param(1, $ditagFeature->hit_strand, SQL_TINYINT);
  $sth->bind_param(1, $ditagFeature->cigar_line, SQL_LONGVARCHAR);
  $sth->bind_param(1, $ditagFeature->ditag_side, SQL_VARCHAR);
  $sth->bind_param(1, $ditagFeature->ditag_pair_id, SQL_INTEGER);
  $sth->bind_param(1, $ditagFeature->dbID, SQL_INTEGER);

  my $result =$sth->execute();

  return $result;
}


=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$dfa->list_dbIDs()};
  Description: Gets an array of internal IDs for all DitagFeature objects in
               the current database.
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : List of ints
  Exceptions : None
  Status     : Stable

=cut

sub list_dbIDs {
  my ($self, $ordered) = @_;
	
  return $self->_list_dbIDs('ditag_feature', undef, $ordered);
}

1;

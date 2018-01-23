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

Bio::EnsEMBL::DBSQL::ExonAdaptor - An adaptor responsible for the retrieval and
storage of exon objects

=head1 SYNOPSIS

  my $exon_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Exon' );

  my $exon = $exon_adaptor->fetch_by_dbID($dbID);

=head1 DESCRIPTION

The ExonAdaptor is responsible for retrieving and storing Exon objects
from an Ensembl database.  Most of the ExonAdaptor functionality is
inherited from the B<Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor> class.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::ExonAdaptor;

use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
 
use vars qw( @ISA );
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

sub _tables {
  my $self = shift;

  # Allow the table definition to be overridden by certain methods.
  if ( defined( $self->{'tables'} ) ) {
    return @{ $self->{'tables'} };
  }

  return ( [ 'exon', 'e' ] );
}


# _columns
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

  my $created_date =
    $self->db->dbc->from_date_to_seconds("created_date");
  my $modified_date =
    $self->db->dbc->from_date_to_seconds("modified_date");

  return (
    'e.exon_id',        'e.seq_region_id',     'e.seq_region_start',
    'e.seq_region_end', 'e.seq_region_strand', 'e.phase',
    'e.end_phase',      'e.is_current',        'e.is_constitutive',
    'e.stable_id',      'e.version',           $created_date,
    $modified_date
  );
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
  my $self = shift;
  return $self->{'final_clause'} || '';
}


sub fetch_all {
  my ($self) = @_;

  my $constraint = 'e.is_current = 1';
  my @exons  = @{ $self->generic_fetch($constraint) };
  return \@exons ;
}

=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id
               the stable id of the exon to retrieve
  Example    : $exon = $exon_adaptor->fetch_by_stable_id('ENSE0000988221');
  Description: Retrieves an Exon from the database via its stable id
  Returntype : Bio::EnsEMBL::Exon in native coordinates.
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
  my ($self, $stable_id) = @_;

  my $constraint = "e.stable_id = ? AND e.is_current = 1";

  $self->bind_param_generic_fetch($stable_id,SQL_VARCHAR);
  my ($exon) = @{ $self->generic_fetch($constraint) };

  # If we didn't get anything back, desperately try to see if there's
  # a version number in the stable_id
  if(!defined($exon) && (my $vindex = rindex($stable_id, '.'))) {
      $exon = $self->fetch_by_stable_id_version(substr($stable_id,0,$vindex),
						substr($stable_id,$vindex+1));
  }

  return $exon;
}


=head2 fetch_by_stable_id_version

  Arg [1]    : String $id 
               The stable ID of the exon to retrieve
  Arg [2]    : Integer $version
               The version of the stable_id to retrieve
  Example    : $exon = $exon_adaptor->fetch_by_stable_id('ENSE0000988221', 3);
  Description: Retrieves an exon object from the database via its stable id and version.
               The exon will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Exon or undef
  Exceptions : if we cant get the exon in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id_version {
    my ($self, $stable_id, $version) = @_;

    # Enforce that version be numeric
    return unless($version =~ /^\d+$/);

    my $constraint = "e.stable_id = ? AND e.version = ? AND e.is_current = 1";
    $self->bind_param_generic_fetch($stable_id, SQL_VARCHAR);
    $self->bind_param_generic_fetch($version, SQL_INTEGER);
    my ($exon) = @{$self->generic_fetch($constraint)};

    return $exon;
}

=head2 fetch_all_versions_by_stable_id 

  Arg [1]     : String $stable_id 
                The stable ID of the exon to retrieve
  Example     : my $exon = $exon_adaptor->fetch_all_version_by_stable_id
                  ('ENSE00000309301');
  Description : Similar to fetch_by_stable_id, but retrieves all versions of an
                exon stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Exon objects
  Exceptions  : if we cant get the gene in given coord system
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_versions_by_stable_id {
  my ($self, $stable_id) = @_;

  my $constraint = "e.stable_id = ?";

  $self->bind_param_generic_fetch($stable_id,SQL_VARCHAR);

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Example    : none
  Description: Retrieves all Exons for the Transcript in 5-3 order
  Returntype : listref Bio::EnsEMBL::Exon on Transcript slice 
  Exceptions : throws if transcript has no slice
  Caller     : Transcript->get_all_Exons()
  Status     : Stable

=cut

sub fetch_all_by_Transcript {
  my ( $self, $transcript ) = @_;

  my $tslice = $transcript->slice();
  my $slice;

  if ( !defined($tslice) ) {
    throw("Transcript must have attached slice to retrieve exons.");
  }

  # use a small slice the same size as the transcript
  if ( !$tslice->is_circular() ) {
    $slice =
      $self->db()->get_SliceAdaptor()->fetch_by_Feature($transcript);
  } else {
    # Circular.
    $slice = $tslice;
  }

  # Override the tables definition to provide an additional join to the
  # exon_transcript table.  For efficiency we cannot afford to have this
  # in as a left join every time.
  my @tables = $self->_tables();

  # Be extra cautious so that we do not add 'exon_transcript' twice.
  my $found = 0;
  foreach my $table (@tables) {
    if ( $table->[0] eq 'exon_transcript' ) {
      $found = 1;
      last;
    }
  }
  if ( !$found ) {
    push @tables, [ 'exon_transcript', 'et' ];
  }

  $self->{'tables'}       = \@tables;
  $self->{'final_clause'} = "ORDER BY et.transcript_id, et.rank";

  my $constraint =
      "et.transcript_id = "
    . $transcript->dbID()
    . " AND e.exon_id = et.exon_id";

  # fetch all of the exons
  my $exons = $self->fetch_all_by_Slice_constraint($slice, $constraint);

  # un-override the table definition
  delete( $self->{'tables'} );
  delete( $self->{'final_clause'} );

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

  Arg [1]    : Bio::EnsEMBL::Exon $exon
               the exon to store in this database
  Example    : $exon_adaptor->store($exon);
  Description: Stores an exon in the database
  Returntype : none
  Exceptions : thrown if exon (or component exons) do not have a contig_id
               or if $exon->start, $exon->end, $exon->strand, or $exon->phase 
               are not defined or if $exon is not a Bio::EnsEMBL::Exon
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ($self, $exon) = @_;

  if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
    throw("$exon is not a EnsEMBL exon - not storing.");
  }

  my $db = $self->db();

  if($exon->is_stored($db)) {
    return $exon->dbID();
  }

  if( ! $exon->start || ! $exon->end ||
      ! $exon->strand || ! defined $exon->phase ) {
    throw("Exon does not have all attributes to store");
  }

  # Default to is_current = 1 if this attribute is not set
  my $is_current = $exon->is_current();
  if ( !defined($is_current) ) { $is_current = 1 }

  # Default to is_constitutive = 0 if this attribute is not set
  my $is_constitutive = $exon->is_constitutive();
  if ( !defined($is_constitutive) ) { $is_constitutive = 0 }

  my $exon_sql = q{
    INSERT into exon ( seq_region_id, seq_region_start,
		       seq_region_end, seq_region_strand, phase,
		       end_phase, is_current, is_constitutive                      
  };
  if ( defined($exon->stable_id) ) {
      my $created = $self->db->dbc->from_seconds_to_date($exon->created_date());
      my $modified = $self->db->dbc->from_seconds_to_date($exon->modified_date());

      $exon_sql .= sprintf ", stable_id, version%s%s) VALUES ( ?,?,?,?,?,?,?,?,?,?%s%s)", $created?", created_date":'', $modified?", modified_date":'', $created?",$created":'', $modified?",$modified":'';
     
  } else {
      $exon_sql .= q{
         ) VALUES ( ?,?,?,?,?,?,?,?)
      };
  }


  my $exonst = $self->prepare($exon_sql);

  my $exonId = undef;

  my $original = $exon;
  my $seq_region_id;
  ($exon, $seq_region_id) = $self->_pre_store($exon);

  #store the exon
  $exonst->bind_param( 1, $seq_region_id,   SQL_INTEGER );
  $exonst->bind_param( 2, $exon->start,     SQL_INTEGER );
  $exonst->bind_param( 3, $exon->end,       SQL_INTEGER );
  $exonst->bind_param( 4, $exon->strand,    SQL_TINYINT );
  $exonst->bind_param( 5, $exon->phase,     SQL_TINYINT );
  $exonst->bind_param( 6, $exon->end_phase, SQL_TINYINT );
  $exonst->bind_param( 7, $is_current,      SQL_TINYINT );
  $exonst->bind_param( 8, $is_constitutive, SQL_TINYINT );

  if ( defined($exon->stable_id) ) {

     $exonst->bind_param(  9, $exon->stable_id, SQL_VARCHAR );
     $exonst->bind_param( 10, $exon->version,   SQL_INTEGER ); 
  }

  $exonst->execute();
  $exonId = $self->last_insert_id('exon_id', undef, 'exon');

  # Now the supporting evidence
  my $esf_adaptor = $db->get_SupportingFeatureAdaptor;
  $esf_adaptor->store($exonId, $exon->get_all_supporting_features);

  #
  # Finally, update the dbID and adaptor of the exon (and any component exons)
  # to point to the new database
  #

  $original->adaptor($self);
  $original->dbID($exonId);

  return $exonId;
}


=head2 remove

  Arg [1]    : Bio::EnsEMBL::Exon $exon
               the exon to remove from the database
  Example    : $exon_adaptor->remove($exon);
  Description: Removes an exon from the database.  This method is generally
               called by the TranscriptAdaptor::store method. Database
               integrity will not be maintained if this method is simply
               called on its own without taking into account transcripts which
               may refer to the exon being removed.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $exon = shift;

  if(!ref($exon) || !$exon->isa('Bio::EnsEMBL::Exon')) {
    throw('Bio::EnsEMBL::Exon argument expected.');
  }

  if(!$exon->is_stored($self->db())) {
    warning("Cannot remove exon " .$exon->dbID.
            "Is not stored in this database.");
    return;
  }

  # sanity check: make sure nobdody tries to slip past a prediction exon
  # which inherits from exon but actually uses different tables
  if($exon->isa('Bio::EnsEMBL::PredictionExon')) {
    throw("ExonAdaptor can only remove Exons not PredictionExons.");
  }

  # Remove the supporting features of this exon

  my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
  my $dna_adp = $self->db->get_DnaAlignFeatureAdaptor;

  my $sth = $self->prepare("SELECT feature_type, feature_id  " .
                           "FROM supporting_feature " .            
			   "WHERE exon_id = ?");
  $sth->bind_param(1, $exon->dbID, SQL_INTEGER);
  $sth->execute();

  # statements to check for shared align_features
  my $sth1 = $self->prepare("SELECT count(*) FROM supporting_feature " .
			    "WHERE feature_type = ? AND feature_id = ?");
  my $sth2 = $self->prepare("SELECT count(*) " .
                            "FROM transcript_supporting_feature " .
			    "WHERE feature_type = ? AND feature_id = ?");

  SUPPORTING_FEATURE:
  while(my ($type, $feature_id) = $sth->fetchrow()){
    
    # only remove align_feature if this is the last reference to it
    $sth1->bind_param(1, $type, SQL_VARCHAR);
    $sth1->bind_param(2, $feature_id, SQL_INTEGER);
    $sth1->execute;
    $sth2->bind_param(1, $type, SQL_VARCHAR);
    $sth2->bind_param(2, $feature_id, SQL_INTEGER);
    $sth2->execute;
    my ($count1) = $sth1->fetchrow;
    my ($count2) = $sth2->fetchrow;
    if ($count1 + $count2 > 1) {
      #warn "shared feature, not removing $type|$feature_id\n";
      next SUPPORTING_FEATURE;
    }
    
    #warn "removing $type|$feature_id\n";
  
    if($type eq 'protein_align_feature'){
      my $f = $prot_adp->fetch_by_dbID($feature_id);
      $prot_adp->remove($f);
    }
    elsif($type eq 'dna_align_feature'){
      my $f = $dna_adp->fetch_by_dbID($feature_id);
      $dna_adp->remove($f);
    }
    else {
      warning("Unknown supporting feature type $type. Not removing feature.");
    }
  }
  $sth->finish();
  $sth1->finish();
  $sth2->finish();

  # delete the association to supporting features

  $sth = $self->prepare("DELETE FROM supporting_feature WHERE exon_id = ?");
  $sth->bind_param(1, $exon->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();


  # delete the exon

  $sth = $self->prepare( "DELETE FROM exon WHERE exon_id = ?" );
  $sth->bind_param(1, $exon->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  $exon->dbID(undef);
  $exon->adaptor(undef);

  return;
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
   my ($self, $ordered) = @_;

   return $self->_list_dbIDs("exon",undef, $ordered);
}


=head2 list_stable_ids

  Arg [1]    : none
  Example    : @stable_exon_ids = @{$exon_adaptor->list_stable_dbIDs()};
  Description: Gets an array of stable ids for all exons in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_stable_ids {
   my ($self) = @_;

   return $self->_list_dbIDs("exon", "stable_id");
}

#_objs_from_sth
#
#  Arg [1]    : StatementHandle $sth
#  Example    : none 
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of Exons
#  Returntype : listref of Bio::EnsEMBL::Exons in target coordinate system
#  Exceptions : none
#  Caller     : internal

sub _objs_from_sth {
  my ( $self, $sth, $mapper, $dest_slice ) = @_;

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

  my ( $exon_id,        $seq_region_id,     $seq_region_start,
       $seq_region_end, $seq_region_strand, $phase,
       $end_phase,      $is_current,        $is_constitutive,
       $stable_id,      $version,           $created_date,
       $modified_date );

  $sth->bind_columns( \( $exon_id, $seq_region_id, $seq_region_start,
                         $seq_region_end, $seq_region_strand, $phase,
                         $end_phase, $is_current, $is_constitutive,
                         $stable_id, $version,    $created_date,
                         $modified_date ) );

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

    # Finally, create the new exon.
    push( @exons,
          $self->_create_feature_fast(
                           'Bio::EnsEMBL::Exon', {
                             'start'         => $seq_region_start,
                             'end'           => $seq_region_end,
                             'strand'        => $seq_region_strand,
                             'adaptor'       => $self,
                             'slice'         => $slice,
                             'dbID'          => $exon_id,
                             'stable_id'     => $stable_id,
                             'version'       => $version,
                             'created_date'  => $created_date || undef,
                             'modified_date' => $modified_date || undef,
                             'phase'         => $phase,
                             'end_phase'     => $end_phase,
                             'is_current'    => $is_current,
                             'is_constitutive' => $is_constitutive } )
    );

  } ## end while ( $sth->fetch() )

  return \@exons;
} ## end sub _objs_from_sth


1;



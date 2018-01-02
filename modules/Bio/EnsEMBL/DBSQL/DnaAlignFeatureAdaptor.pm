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

Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor - Adaptor for DnaAlignFeatures

=head1 SYNOPSIS

  $dafa = $registry->get_adaptor( 'Human', 'Core', 'DnaAlignFeature' );

  my @features = @{ $dafa->fetch_all_by_Slice($slice) };

  $dafa->store(@features);

=head1 DESCRIPTION

This is an adaptor responsible for the retrieval and storage of
DnaDnaAlignFeatures from the database. This adaptor inherits most of its
functionality from the BaseAlignFeatureAdaptor and BaseFeatureAdaptor
superclasses.

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor);


=head2 _tables

  Args       : none
  Example    : @tabs = $self->_tables
  Description: PROTECTED implementation of the abstract method inherited from
               BaseFeatureAdaptor.  Returns list of [tablename, alias] pairs
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::generic_fetch
  Status     : Stable

=cut

sub _tables {
  my $self = shift;

  return (['dna_align_feature', 'daf'],['external_db','exdb']);
}


sub _left_join{
    return (['external_db',"exdb.external_db_id = daf.external_db_id"]);
}

=head2 _columns

  Args       : none
  Example    : @columns = $self->_columns
  Description: PROTECTED implementation of abstract superclass method.  
               Returns a list of columns that are needed for object creation.
  Returntype : list of strings
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::generic_fetch
  Status     : Stable

=cut

sub _columns {
  my $self = shift;

  #warning, implementation of _objs_from_sth method depends on order of list
  return qw(daf.dna_align_feature_id
            daf.seq_region_id
            daf.analysis_id
            daf.seq_region_start
            daf.seq_region_end
            daf.seq_region_strand
            daf.hit_start
            daf.hit_end
            daf.hit_name
            daf.hit_strand
            daf.cigar_line
            daf.evalue
            daf.perc_ident
            daf.score
            daf.external_db_id
            daf.hcoverage
            daf.align_type
            exdb.db_name
            exdb.db_display_name);
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::DnaAlignFeatures @feats
               the features to store in the database
  Example    : $dna_align_feature_adaptor->store(@features);
  Description: Stores a list of DnaAlignFeatures in the database
  Returntype : none
  Exceptions : throw if any of the provided features cannot be stored
               which may occur if:
                 * The feature does not have an associate Slice
                 * The feature does not have an associated analysis
                 * The Slice the feature is associated with is on a seq_region
                   unknown to this database
               A warning is given if:
                 * The feature has already been stored in this db
  Caller     : Pipeline
  Status     : Stable

=cut

sub store {
  my ( $self, @feats ) = @_;

  throw("Must call store with features") if ( scalar(@feats) == 0 );

  my @tabs = $self->_tables;
  my ($tablename) = @{ $tabs[0] };

  my $db               = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

  my $sth = $self->prepare(
    "INSERT INTO $tablename (seq_region_id, seq_region_start,
                             seq_region_end, seq_region_strand,
                             hit_start, hit_end, hit_strand, hit_name,
                             cigar_line, analysis_id, score, evalue,
                             perc_ident, external_db_id, hcoverage, align_type)
     VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"    # 16 arguments, external_data removed and align_type added
  );

FEATURE:
  foreach my $feat (@feats) {
    if ( !ref $feat || !$feat->isa("Bio::EnsEMBL::DnaDnaAlignFeature") )
    {
      throw("feature must be a Bio::EnsEMBL::DnaDnaAlignFeature, not a [". ref($feat). "]." );
    }

    if ( $feat->is_stored($db) ) {
      warning( "DnaDnaAlignFeature [".$feat->dbID()."] is already stored in this database." );
      next FEATURE;
    }

    my $hstart  = $feat->hstart();
    my $hend    = $feat->hend();
    my $hstrand = $feat->hstrand();
    $self->_check_start_end_strand( $hstart, $hend, $hstrand, $feat->hslice );

    my $cigar_string = $feat->cigar_string();
    if ( !$cigar_string ) {
      $cigar_string = $feat->length() . 'M';
      warning( "DnaDnaAlignFeature does not define a cigar_string.\n"
          . "Assuming ungapped block with cigar_line=$cigar_string ." );
    }
    my $align_type = $feat->align_type();

    my $hseqname = $feat->hseqname();
    if ( !$hseqname ) {
      throw("DnaDnaAlignFeature must define an hseqname.");
    }

    if ( !defined( $feat->analysis ) ) {
      throw(
        "An analysis must be attached to the features to be stored.");
    }

    #store the analysis if it has not been stored yet
    if ( !$feat->analysis->is_stored($db) ) {
      $analysis_adaptor->store( $feat->analysis() );
    }

    my $original = $feat;
    my $seq_region_id;
    ( $feat, $seq_region_id ) = $self->_pre_store($feat);
    
    $sth->bind_param( 1,  $seq_region_id,        SQL_INTEGER );
    $sth->bind_param( 2,  $feat->start,          SQL_INTEGER );
    $sth->bind_param( 3,  $feat->end,            SQL_INTEGER );
    $sth->bind_param( 4,  $feat->strand,         SQL_TINYINT );
    $sth->bind_param( 5,  $hstart,               SQL_INTEGER );
    $sth->bind_param( 6,  $hend,                 SQL_INTEGER );
    $sth->bind_param( 7,  $hstrand,              SQL_TINYINT );
    $sth->bind_param( 8,  $hseqname,             SQL_VARCHAR );
    $sth->bind_param( 9,  $cigar_string,         SQL_LONGVARCHAR );
    $sth->bind_param( 10, $feat->analysis->dbID, SQL_INTEGER );
    $sth->bind_param( 11, $feat->score,          SQL_DOUBLE );
    $sth->bind_param( 12, $feat->p_value,        SQL_DOUBLE );
    $sth->bind_param( 13, $feat->percent_id,     SQL_FLOAT );
    $sth->bind_param( 14, $feat->external_db_id, SQL_INTEGER );
    $sth->bind_param( 15, $feat->hcoverage,      SQL_DOUBLE );
    $sth->bind_param( 16, $feat->align_type,     SQL_VARCHAR);
    $sth->execute();

    my $dbId = $self->last_insert_id("${tablename}_id", undef, $tablename);

    # store attributes if there are any
    my $attr_adaptor = $db->get_AttributeAdaptor();
    $attr_adaptor->store_on_DnaDnaAlignFeature($dbId, $feat->get_all_Attributes);

    $original->dbID( $dbId );
    $original->adaptor($self);
  } ## end foreach my $feat (@feats)

  $sth->finish();
} ## end sub store


=head2 remove

  Arg [1]    : A feature $feature 
  Example    : $feature_adaptor->remove($feature);
  Description: This removes a feature from the database.  The table the
               feature is removed from is defined by the abstract method
               _tablename, and the primary key of the table is assumed
               to be _tablename() . '_id'.  The feature argument must 
               be an object implementing the dbID method, and for the
               feature to be removed from the database a dbID value must
               be returned.
  Returntype : none
  Exceptions : thrown if $feature arg does not implement dbID(), or if
               $feature->dbID is not a true value
  Caller     : general
  Status     : Stable

=cut


sub remove {
  my ($self, $feature) = @_;

  if(!$feature || !ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Feature argument is required');
  }

  if(!$feature->is_stored($self->db)) {
    throw("This feature is not stored in this database");
  }

  my @tabs = $self->_tables;
  my ($table) = @{$tabs[0]};

  my $sth = $self->prepare("DELETE FROM $table WHERE ${table}_id = ?");
  $sth->bind_param(1,$feature->dbID,SQL_INTEGER);
  $sth->execute();

  # remove the attributes associated with this feature
  my $attrib_adaptor = $self->db->get_AttributeAdaptor;
  $attrib_adaptor->remove_from_DnaDnaAlignFeature($feature);

  #unset the feature dbID ad adaptor
  $feature->dbID(undef);
  $feature->adaptor(undef);

  return;
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle $sth
               an exectuted DBI statement handle generated by selecting 
               the columns specified by _columns() from the table specified 
               by _table()
  Example    : @dna_dna_align_feats = $self->_obj_from_hashref
  Description: PROTECTED implementation of superclass abstract method. 
               Creates DnaDnaAlignFeature objects from a DBI hashref
  Returntype : listref of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseFeatureAdaptor::generic_fetch
  Status     : Stable

=cut

sub _objs_from_sth {
  my ( $self, $sth, $mapper, $dest_slice ) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  # In case of userdata we need the features on the dest_slice.  In case
  # of get_all_supporting_features dest_slice is not provided.
  my $sa = (   $dest_slice
             ? $dest_slice->adaptor()
             : $self->db()->get_SliceAdaptor() );
  my $aa = $self->db->get_AnalysisAdaptor();

  my @features;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my ( $dna_align_feature_id, $seq_region_id,
       $analysis_id,          $seq_region_start,
       $seq_region_end,       $seq_region_strand,
       $hit_start,            $hit_end,
       $hit_name,             $hit_strand,
       $cigar_line,           $evalue,
       $perc_ident,           $score,
       $external_db_id,       $hcoverage,
       $align_type,
       $external_db_name,     $external_display_db_name );

  $sth->bind_columns( \( $dna_align_feature_id, $seq_region_id,
                         $analysis_id,          $seq_region_start,
                         $seq_region_end,       $seq_region_strand,
                         $hit_start,            $hit_end,
                         $hit_name,             $hit_strand,
                         $cigar_line,           $evalue,
                         $perc_ident,           $score,
                         $external_db_id,       $hcoverage,
                         $align_type,
                         $external_db_name, $external_display_db_name )
  );

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

    #get the analysis object
    my $analysis = $analysis_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);
    $analysis_hash{$analysis_id} = $analysis;

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

    # Finally, create the new DnaAlignFeature.
    push( @features,
          $self->_create_feature_fast(
             'Bio::EnsEMBL::DnaDnaAlignFeature', {
               'slice'           => $slice,
               'start'           => $seq_region_start,
               'end'             => $seq_region_end,
               'strand'          => $seq_region_strand,
               'hseqname'        => $hit_name,
               'hstart'          => $hit_start,
               'hend'            => $hit_end,
               'hstrand'         => $hit_strand,
               'score'           => $score,
               'p_value'         => $evalue,
               'percent_id'      => $perc_ident,
               'cigar_string'    => $cigar_line,
               'analysis'        => $analysis,
               'adaptor'         => $self,
               'dbID'            => $dna_align_feature_id,
               'external_db_id'  => $external_db_id,
               'hcoverage'       => $hcoverage,
               'align_type'      => $align_type,
               'dbname'          => $external_db_name,
               'db_display_name' => $external_display_db_name
             } ) );

  } ## end while ( $sth->fetch() )

  return \@features;
} ## end sub _objs_from_sth

=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$dna_align_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all dna align features in 
               the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self, $ordered) = @_;

   return $self->_list_dbIDs("dna_align_feature",undef, $ordered);
}

1;



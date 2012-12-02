=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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
	    daf.external_data
	    daf.pair_dna_align_feature_id
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
                             perc_ident, external_db_id, hcoverage,
                             pair_dna_align_feature_id)
     VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"    # 16 arguments
  );

FEATURE:
  foreach my $feat (@feats) {
    if ( !ref $feat || !$feat->isa("Bio::EnsEMBL::DnaDnaAlignFeature") )
    {
      throw("feature must be a Bio::EnsEMBL::DnaDnaAlignFeature,"
          . " not a ["
          . ref($feat)
          . "]." );
    }

    if ( $feat->is_stored($db) ) {
      warning( "DnaDnaAlignFeature ["
          . $feat->dbID()
          . "] is already stored in this database." );
      next FEATURE;
    }

    my $hstart  = $feat->hstart();
    my $hend    = $feat->hend();
    my $hstrand = $feat->hstrand();
    $self->_check_start_end_strand( $hstart, $hend, $hstrand );

    my $cigar_string = $feat->cigar_string();
    if ( !$cigar_string ) {
      $cigar_string = $feat->length() . 'M';
      warning( "DnaDnaAlignFeature does not define a cigar_string.\n"
          . "Assuming ungapped block with cigar_line=$cigar_string ." );
    }

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
    $sth->bind_param( 16, $feat->pair_dna_align_feature_id,
      SQL_INTEGER );

    $sth->execute();

    $original->dbID( $sth->{'mysql_insertid'} );
    $original->adaptor($self);
  } ## end foreach my $feat (@feats)

  $sth->finish();
} ## end sub store


sub save {
  my ($self, $features) = @_;

  my @feats = @$features;
  throw("Must call store with features") if( scalar(@feats) == 0 );

  my @tabs = $self->_tables;
  my ($tablename) = @{$tabs[0]};

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

  my $sql = qq{INSERT INTO $tablename (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_start, hit_end, hit_strand, hit_name, cigar_line, analysis_id, score, evalue, perc_ident, external_db_id, hcoverage, pair_dna_align_feature_id, external_data) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)};

  my %analyses = ();

  my $sth = $self->prepare($sql);
     
 FEATURE: foreach my $feat ( @feats ) {
    if( !ref $feat || !$feat->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {
      throw("feature must be a Bio::EnsEMBL::DnaDnaAlignFeature,"
            . " not a [".ref($feat)."].");
    }

    if($feat->is_stored($db)) {
      warning("DnaDnaAlignFeature [".$feat->dbID."] is already stored" .
              " in this database.");
      next FEATURE;
    }

    my $hstart  = $feat->hstart || 0; # defined $feat->hstart  ? $feat->hstart : $feat->start ;
    my $hend    = $feat->hend   || 0; # defined $feat->hend    ? $feat->hend : $feat->end;
    my $hstrand = $feat->hstrand|| 0; # defined $feat->hstrand ? $feat->hstrand : $feat->strand;
    if( $hstart && $hend ) {
      if($hend < $hstart) {
        throw("Invalid Feature start/end [$hstart/$hend]. Start must be less than or equal to end.");
      }
    }
    my $cigar_string = $feat->cigar_string();
    if(!$cigar_string) {
      $cigar_string = $feat->length() . 'M';
      warning("DnaDnaAlignFeature does not define a cigar_string.\n" .
              "Assuming ungapped block with cigar_line=$cigar_string .");
    }

    my $hseqname = $feat->hseqname();
    if(!$hseqname) {
      throw("DnaDnaAlignFeature must define an hseqname.");
    }

    if(!defined($feat->analysis)) {
      throw("An analysis must be attached to the features to be stored.");
    }

    #store the analysis if it has not been stored yet
    if(!$feat->analysis->is_stored($db)) {
      $analysis_adaptor->store($feat->analysis());
    }

    $analyses{ $feat->analysis->dbID }++;

    my $original = $feat;
    my $seq_region_id;
    ($feat, $seq_region_id) = $self->_pre_store_userdata($feat);

    my $extra_data = $feat->extra_data ? $self->dump_data($feat->extra_data) : '';

    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$feat->start,SQL_INTEGER);
    $sth->bind_param(3,$feat->end,SQL_INTEGER);
    $sth->bind_param(4,$feat->strand,SQL_TINYINT);
    $sth->bind_param(5,$hstart,SQL_INTEGER);
    $sth->bind_param(6,$hend,SQL_INTEGER);
    $sth->bind_param(7,$hstrand,SQL_TINYINT);
    $sth->bind_param(8,$hseqname,SQL_VARCHAR);
    $sth->bind_param(9,$cigar_string,SQL_LONGVARCHAR);
    $sth->bind_param(10,$feat->analysis->dbID,SQL_INTEGER);
    $sth->bind_param(11,$feat->score,SQL_DOUBLE);
#    $sth->bind_param(11,$feat->score); # if the above statement does not work it means you need to upgrade DBD::mysql, meantime you can replace it with this line
    $sth->bind_param(12,$feat->p_value,SQL_DOUBLE);
    $sth->bind_param(13,$feat->percent_id,SQL_FLOAT);
    $sth->bind_param(14,$feat->external_db_id,SQL_INTEGER);
    $sth->bind_param(15,$feat->hcoverage,SQL_DOUBLE);
    $sth->bind_param(16,$feat->pair_dna_align_feature_id,SQL_INTEGER);
    $sth->bind_param(17,$extra_data,SQL_LONGVARCHAR);


    $sth->execute();
    $original->dbID($sth->{'mysql_insertid'});
    $original->adaptor($self);
  }

  $sth->finish();

## js5 hack to update meta_coord table... 
  if( keys %analyses ) {

    my $sth = $self->prepare( 'select sr.coord_system_id, max(daf.seq_region_end-daf.seq_region_start) from seq_region as sr, dna_align_feature as daf where daf.seq_region_id=sr.seq_region_id and analysis_id in ('.join(',',keys %analyses).') group by coord_system_id' );
    $sth->execute;

    foreach( @{ $sth->fetchall_arrayref } ) {
      my $sth2 = $self->prepare( qq(insert ignore into meta_coord values("dna_align_feature",$_->[0],$_->[1])) );
      $sth2->execute;
      $sth2->finish;

      $sth2 = $self->prepare( qq(update meta_coord set max_length = $_->[1] where coord_system_id = $_->[0] and table_name="dna_align_feature" and max_length < $_->[1]) );
      $sth2->execute;
      $sth2->finish;
    }

    $sth->finish;
  }

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
       $extra_data,           $pair_dna_align_feature_id,
       $external_db_name,     $external_display_db_name );

  $sth->bind_columns( \( $dna_align_feature_id, $seq_region_id,
                         $analysis_id,          $seq_region_start,
                         $seq_region_end,       $seq_region_strand,
                         $hit_start,            $hit_end,
                         $hit_name,             $hit_strand,
                         $cigar_line,           $evalue,
                         $perc_ident,           $score,
                         $external_db_id,       $hcoverage,
                         $extra_data,       $pair_dna_align_feature_id,
                         $external_db_name, $external_display_db_name )
  );

  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;

  if ( defined($mapper) ) {
    $asm_cs      = $mapper->assembled_CoordSystem();
    $cmp_cs      = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $cmp_cs_vers = $cmp_cs->version();
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  my $dest_slice_sr_name;
  my $dest_slice_seq_region_id;

  if ( defined($dest_slice) ) {
    $dest_slice_start         = $dest_slice->start();
    $dest_slice_end           = $dest_slice->end();
    $dest_slice_strand        = $dest_slice->strand();
    $dest_slice_length        = $dest_slice->length();
    $dest_slice_sr_name       = $dest_slice->seq_region_name();
    $dest_slice_seq_region_id = $dest_slice->get_seq_region_id();
  }

FEATURE:
  while ( $sth->fetch() ) {
    # Get the analysis object.
    my $analysis = $analysis_hash{$analysis_id} ||=
      $aa->fetch_by_dbID($analysis_id);

    # Get the slice object.
    my $slice = $slice_hash{ "ID:" . $seq_region_id };

    if ( !defined($slice) ) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      if ( defined($slice) ) {
        $slice_hash{ "ID:" . $seq_region_id } = $slice;
        $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
        $sr_cs_hash{$seq_region_id}   = $slice->coord_system();
      }
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};

    # Remap the feature coordinates to another coord system
    # if a mapper was provided.
    if ( defined($mapper) ) {

	if (defined $dest_slice && $mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper')  ) {
	    ( $seq_region_id,  $seq_region_start,
	      $seq_region_end, $seq_region_strand )
		=
		$mapper->map( $sr_name, $seq_region_start, $seq_region_end,
                          $seq_region_strand, $sr_cs, 1, $dest_slice);

	} else {

	    ( $seq_region_id,  $seq_region_start,
	      $seq_region_end, $seq_region_strand )
		=
		$mapper->fastmap( $sr_name, $seq_region_start, $seq_region_end,
                          $seq_region_strand, $sr_cs );
	}

      # Skip features that map to gaps or coord system boundaries.
      if ( !defined($seq_region_id) ) { next FEATURE }

      # Get a slice in the coord system we just mapped to.
      if ( $asm_cs == $sr_cs
           || ( $cmp_cs != $sr_cs && $asm_cs->equals($sr_cs) ) )
      {
        $slice = $slice_hash{ "ID:" . $seq_region_id } ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
      } else {
        $slice = $slice_hash{ "ID:" . $seq_region_id } ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
      }
    }

    # If a destination slice was provided, convert the coords.  If the
    # dest_slice starts at 1 and is forward strand, nothing needs doing.
    if ( defined($dest_slice) ) {
      if ( $dest_slice_start != 1 || $dest_slice_strand != 1 ) {
        if ( $dest_slice_strand == 1 ) {
          $seq_region_start = $seq_region_start - $dest_slice_start + 1;
          $seq_region_end   = $seq_region_end - $dest_slice_start + 1;
        } else {
          my $tmp_seq_region_start = $seq_region_start;
          $seq_region_start = $dest_slice_end - $seq_region_end + 1;
          $seq_region_end = $dest_slice_end - $tmp_seq_region_start + 1;
          $seq_region_strand = -$seq_region_strand;
        }

        # Throw away features off the end of the requested slice.
        if (    $seq_region_end < 1
             || $seq_region_start > $dest_slice_length
             || ( $dest_slice_seq_region_id ne $seq_region_id ) )
        {
          next FEATURE;
        }
      }
      $slice = $dest_slice;
    }

    # Inlining the following in the hash causes major issues with 5.16 and messes up the hash 
    my $evalled_extra_data = $extra_data ? $self->get_dumped_data($extra_data) : '';

    # Finally, create the new DnaAlignFeature.
    push( @features,
          $self->_create_feature_fast(
             'Bio::EnsEMBL::DnaDnaAlignFeature', {
               'slice'          => $slice,
               'start'          => $seq_region_start,
               'end'            => $seq_region_end,
               'strand'         => $seq_region_strand,
               'hseqname'       => $hit_name,
               'hstart'         => $hit_start,
               'hend'           => $hit_end,
               'hstrand'        => $hit_strand,
               'score'          => $score,
               'p_value'        => $evalue,
               'percent_id'     => $perc_ident,
               'cigar_string'   => $cigar_line,
               'analysis'       => $analysis,
               'adaptor'        => $self,
               'dbID'           => $dna_align_feature_id,
               'external_db_id' => $external_db_id,
               'hcoverage'      => $hcoverage,
               'extra_data'     => $evalled_extra_data,
               'dbname'                    => $external_db_name,
               'db_display_name'           => $external_display_db_name,
               'pair_dna_align_feature_id' => $pair_dna_align_feature_id
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



#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::AffyFeatureAdaptor
#
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AffyFeatureAdaptor

=head1 SYNOPSIS

  my $featureAdaptor = $db->get_AffyFeatureAdaptor();
  
  my $features = $featureAdaptor->fetch_all_by_AffyProbe( $probe );
  $features = $featureAdaptor->fetch_all_by_Slice_arrayname( $slice, "Array1", "Array2" );

  
=head1 DESCRIPTION
 
  Adaptor to load and store AffyFeature objects in the connected database. 
  The store does not look for duplicates, so you have to make sure to store 
  Features only once.
 
=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::AffyFeatureAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::AffyFeature;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 fetch_all_by_AffyProbe

  Arg [1]    : Bio::EnsEMBL::AffyProbe $probe
  Example    : none
  Description: Fetches all featurs that given probe creates.
  Returntype : listref of Bio::EnsEMBL::AffyFeature
  Exceptions : none
  Caller     : AffyProbe->get_all_AffyFeatures()
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut


sub fetch_all_by_AffyProbe {
  my $self = shift;
  my $probe = shift;

  throw( "AffyProbe is required argument to fetch_by_AffyProbe\nYou supplied $probe" ) unless 
    ref($probe) && $probe->isa( "Bio::EnsEMBL::AffyProbe" );

  throw "Supplied probe is not stored in database" unless defined $probe->dbID();

  return $self->generic_fetch( "af.affy_probe_id = ". $probe->dbID() );
}


=head2 fetch_all_by_Slice_arrayname

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2..]  : list of String @arraynames
  Example    : my $slice = $sliceAdaptor->fetch_by_region( "chromosome", "1" );
               $features = $affy_feature_adaptor->fetch_by_Slice_arrayname( $slice, "" );
  Description: Retrieves a list of features on given Slice that are created by probes from the
               Arrays that are given by their names. If none is given, all features on that Slice 
               are returned.
  Returntype : listref of Bio::EnsEMBL::AffyFeature
  Exceptions : throws if arraynames passed
  Caller     : slice->get_all_AffyFeatures()
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub fetch_all_by_Slice_arrayname {
  my( $self, $slice, @arraynames ) = @_;

  throw "Need arrayname as parameter" unless @arraynames;

  my $constraint;

  if( scalar @arraynames == 1 ) {
    $constraint = qq( aa.name = "$arraynames[0]" );
  } else { 
    $constraint = join '","' ,  @arraynames;
    $constraint = qq( aa.name in ( "$constraint" ) );
  }

  return $self->SUPER::fetch_all_by_Slice_constraint( $slice, $constraint );
}

=head2 store

  Arg [1]    : list of Bio::EnsEMBL::AffyFeatures @afs
               the affy features to store in the database
  Example    : $affy_feature_adaptor->store(@affy_feats);
  Description: Stores a list of affy feature objects in the database
  Returntype : none
  Exceptions : thrown if @afs is not defined, if any of the features do not
               have an attached slice or atached probe.
               or if any elements of @afs are not Bio::EnsEMBL::AffyFeatures 
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub store{
  my ($self,@afs) = @_;

  if( scalar(@afs) == 0 ) {
    throw("Must call store with list of AffyFeatures");
  }

  my $sth = $self->prepare
    ("INSERT INTO affy_feature (seq_region_id, seq_region_start, " .
     "                          seq_region_end, seq_region_strand, " .
     "                          affy_probe_id, analysis_id, mismatches ) " .
     "VALUES (?,?,?,?,?,?,?)" );

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

 FEATURE: foreach my $af ( @afs ) {

    if( !ref $af || !$af->isa("Bio::EnsEMBL::AffyFeature") ) {
      throw("Feature must be an Ensembl AffyFeature, " .
            "not a [".ref($af)."]");
    }

    if($af->is_stored($db)) {
      warning("AffyFeature [".$af->dbID."] is already stored" .
              " in this database.");
      next FEATURE;
    }

    if(!defined($af->analysis)) {
      throw("An analysis must be attached to the features to be stored.");
    }

    #store the analysis if it has not been stored yet
    if(!$af->analysis->is_stored($db)) {
      $analysis_adaptor->store($af->analysis());
    }

    my $original = $af;
    my $seq_region_id;
    ($af, $seq_region_id) = $self->_pre_store($af);

    $sth->execute($seq_region_id, $af->start, $af->end, $af->strand,
                  $af->probe->dbID, $af->analysis->dbID, $af->mismatchcount );

    $original->dbID($sth->{'mysql_insertid'});
    $original->adaptor($self);
  }
}


=head2 _tablename

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub _tables {
  my $self = shift;
  
  return ( [ 'affy_feature', 'af' ], 
	   [ 'affy_probe',   'ap' ], 
	   [ 'affy_array',   'aa' ]);
}


sub _default_where_clause {
  my $self = shift;

  return 'af.affy_probe_id = ap.affy_probe_id and ap.affy_array_id=aa.affy_array_id';
}

=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a list of columns to use for queries
  Returntype : list of strings
  Exceptions : none
  Caller     : internal
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub _columns {
  my $self = shift;

  return qw( af.affy_feature_id af.seq_region_id af.seq_region_start
             af.seq_region_end af.seq_region_strand 
             af.mismatches af.affy_probe_id af.analysis_id
             aa.name ap.probeset ap.name)

}


=head2 _objs_from_sth

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: PROTECTED implementation of superclass abstract method.
               creates SimpleFeatures from an executed DBI statement handle.
  Returntype : list reference to Bio::EnsEMBL::AffyFeature objects
  Exceptions : none
  Caller     : internal
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.

  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db->get_AnalysisAdaptor();

  my @features;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;


  my( $affy_feature_id,$seq_region_id, $seq_region_start, $seq_region_end,
      $seq_region_strand, $mismatches, $analysis_id,
      $affy_probe_id, $array_name, $affy_probe_set, $affy_probe_name );

  $sth->bind_columns(\$affy_feature_id,\$seq_region_id, \$seq_region_start,
                     \$seq_region_end, \$seq_region_strand, 
		     \$mismatches, \$affy_probe_id, \$analysis_id,
                     \$array_name, \$affy_probe_set, \$affy_probe_name );
 
  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;
  if($mapper) {
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
  if( $dest_slice ) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
  }

  my $last_feature_id = -1;
  FEATURE: while($sth->fetch()) {
    #get the analysis object
    next if( $last_feature_id == $affy_feature_id );
    $last_feature_id = $affy_feature_id;

    my $analysis = $analysis_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);

    #get the slice object
    my $slice = $slice_hash{"ID:".$seq_region_id};

    unless($slice) {
      $slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id}       = $slice->coord_system();
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};
 
   #
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    #
    if($mapper) {

      ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
        $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($sr_name));

      #get a slice in the coord system we just mapped to
      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"} ||=
          $sa->fetch_by_region($cmp_cs_name, $sr_name,undef, undef, undef, $cmp_cs_vers);
      } else {
        $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
          $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef, $asm_cs_vers);
      }
    }

    #
    # If a destination slice was provided convert the coords
    # If the dest_slice starts at 1 and is foward strand, nothing needs doing
    #
    if($dest_slice) {
      unless( $dest_slice_start == 1 && $dest_slice_strand == 1) {
        if($dest_slice_strand == 1) {
          $seq_region_start = $seq_region_start - $dest_slice_start + 1;
          $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
        } else {
          my $tmp_seq_region_start = $seq_region_start;
          $seq_region_start        = $dest_slice_end - $seq_region_end       + 1;
          $seq_region_end          = $dest_slice_end - $tmp_seq_region_start + 1;
          $seq_region_strand      *= -1;
        }
      }

      #throw away features off the end of the requested slice
      next FEATURE if $seq_region_end < 1 || $seq_region_start > $dest_slice_length ||
	  ( $dest_slice_sr_name ne $sr_name );

      $slice = $dest_slice;
    }

    push @features, Bio::EnsEMBL::AffyFeature->new_fast( {
      'start'         => $seq_region_start,
      'end'           => $seq_region_end,
      'strand'        => $seq_region_strand,
      'slice'         => $slice,
      'analysis'      => $analysis,
      'adaptor'       => $self,
      'dbID'          => $affy_feature_id,
      'mismatchcount' => $mismatches,
      '_probe_id'     => $affy_probe_id,
      'probeset'      => $affy_probe_set,
      '_probe_name'   => $affy_probe_name
    } );
  }

  return \@features;
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$simple_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub list_dbIDs {
  my ($self) = @_;
  return $self->_list_dbIDs( "affy_feature" );
}

1;

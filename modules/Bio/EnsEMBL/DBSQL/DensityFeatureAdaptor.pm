#
# Ensembl module for Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor

=head1 SYNOPSIS

my $dfa = $database_adaptor->get_DensityFeatureAdaptor();

my $interpolate = 1;
my $blocks_wanted = 50;

@dense_feats = @{$dfa->fetch_all_by_Slice($slice,'SNPDensity',
                                             $blocks_wanted, $interpolate);}

=head1 DESCRIPTION

Density Feature Adaptor - An adaptor responsible for the creation of density
features from the database.

=head1 CONTACT

Post questions to the Ensembl developer list.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor;
use vars qw(@ISA);
use strict;


use POSIX;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);



=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice - The slice representing the region
               to retrieve density features from.
  Arg [2]    : string $logic_name - The logic name of the density features to
               retrieve.
  Arg [3]    : int $num_blocks (optional; default = 50) - The number of
               features that are desired. The ratio between the size of these
               features and the size of the features in the database will be
               used to determine which database features will be used.
  Arg [4]    : boolean $interpolate (optional; default = 0) - A flag indicating
               whether the features in the database should be interpolated to
               fit them to the requested number of features.  If true the
               features will be interpolated to provide $num_blocks features.
               This will not guarantee that exactly $num_blocks features are
               returned due to rounding etc. but it will be close.
  Arg [5]    : float $max_ratio - The maximum ratio between the size of the
               requested features (as determined by $num_blocks) and the actual
               size of the features in the database.  If this value is exceeded
               then an empty list will be returned.  This can be used to
               prevent excessive interpolation of the database values.
  Example    : #interpolate:
               $feats = $dfa->fetch_all_by_Slice($slice,'SNPDensity', 10, 1);
               #do not interpoloate, get what is in the database:
               $feats = $dfa->fetch_all_by_Slice($slice,'SNPDensity', 50);
               #interpolate, but not too much
               $feats = $dfa->fetch_all_by_Slice($slice,'SNPDensity',50,1,5.0);
  Description: Retrieves a set of density features which overlap the region
               of this slice. Density features are a discrete representation
               of a continuous value along a sequence, such as a density or
               percent coverage.  Density Features may be stored in chunks of
               different sizes in the database, and interpolated to fit the
               desired size for the requested region.  For example the database
               may store a single density value for each 1KB and also for each
               1MB.  When fetching for an entire chromosome the 1MB density
               chunks will be used if the requested number of blocks is not
               very high.
  Returntype : Bio::EnsEMBL::DensityFeature
  Exceptions : warning on invalid $num_blocks argument
               warning if no features with logic_name $logic_name exist
               warning if density_type table has invalid block_size value
  Caller     : general

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $logic_name, $num_blocks, $interpolate, $max_ratio) = @_;

  if(defined($num_blocks) && $num_blocks < 1) {
    warning("Invalid number of density blocks [$num_blocks] requested.\n" .
	   "Returning empty list.");
    return [];
  }

  $num_blocks ||= 50;
  my $length = $slice->length();


#  my $wanted_block_size = int($length/$num_blocks);
  my $wanted_block_size = POSIX::ceil($length/$num_blocks);

  my $analysis_adaptor = $self->db()->get_AnalysisAdaptor();
  my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

  if(!$analysis) {
    warning("No Analysis with logic_name=[$logic_name] found.\n" .
            'Returning empty list.');
    return [];
  }

  my $sth = $self->prepare
    ("SELECT density_type_id, block_size, analysis_id, value_type  ".
     "FROM   density_type " .
     "WHERE  analysis_id = ?");
  $sth->execute($analysis->dbID());

  my $best_ratio = undef;
  my ($density_type_id, $density_value_type, $new_block_size, $row);

  while($row = $sth->fetchrow_arrayref) {
    my $block_size = $row->[1];

    if($block_size < 1) {
      warning("density_type table contains invalid block_size=$block_size.");
      next;
    }

    my $ratio;
    if($block_size < $wanted_block_size) {
      $ratio = $wanted_block_size/$block_size;
    } else {
      $ratio = $block_size/$wanted_block_size;
    }

    if(!defined($best_ratio) || $ratio < $best_ratio) {
      $best_ratio = $ratio;
      $density_type_id = $row->[0];
      $new_block_size = $row->[1];
      $density_value_type = $row->[3];
    }
  }

  #the ratio was not good enough, or this logic name was not in the
  #density type table, return an empty list
  if(!defined($best_ratio) ||
     (defined($max_ratio) && $best_ratio > $max_ratio)) {
    return [];
  }


  my $density_type = Bio::EnsEMBL::DensityType->new(-analysis => $analysis,
						    -block_size =>  $new_block_size,
						    -value_type => $density_value_type);

  $density_type->dbID($density_type_id);

  my @features =
    @{$self->fetch_all_by_slice_and_density_type($slice,$density_type)};

  #we don't want to interpolate if the ratio was very close
  $interpolate = 0 if($best_ratio < 1.05);

  return \@features if(!$interpolate);

  #interpolate the features into new features of a different size
  my @out;
  #sort the features on start position
  @features = sort({$a->start() <=> $b->end()} @features);

  #resize the features that were returned
  my $start = 1;
  my $end   = $start+$wanted_block_size-1;

  while($start < $length) {
#    $end = $length if($end > $length);

    my $density_value = 0.0;
    my ($f, $fstart, $fend, $portion);
    my @dvalues;

    #construct a new feature using all of the old density features that
    #we overlapped
    while(($f = shift(@features)) && $end > $f->{'start'}) {

      #what portion of this feature are we using to construct our new block?
      $fstart = ($f->{'start'} < $start) ? $start : $f->{'start'};
      $fend   = ($f->{'end'}   > $end  ) ? $end   : $f->{'end'};
      $fend   = $length if($fend > $length);

      if($density_value_type eq 'sum') {

        $portion = ($fend-$fstart+1)/$f->length();

        #take a percentage of density value, depending on how much of the
        #feature we overlapped
        $density_value += $portion * $f->{'density_value'};

      } elsif($density_value_type eq 'ratio') {

        #maintain a running total of the length used from each feature
        #and its value
        push(@dvalues,[$fend-$fstart+1,$f->{'density_value'}]);

      } else {
        throw("Unknown density value type [$density_value_type].");
      }

      #do not want to look at next feature if we only used part of this one:
      last if($fend < $f->{'end'});
    }

    #if we did not completely overlap the last feature, put it back on so
    #it can be partially used by the next block
    if(defined($f) && (!defined($fend) || $fend < $f->{'end'})) {
      unshift(@features, $f);
    }

    if($density_value_type eq 'ratio') {
      #take a weighted average of the all the density values of the features
      #used to construct this one
      my $total_len = $end - $start + 1;
      if($total_len > 0) {
        foreach my $pair (@dvalues) {
          my ($dlen, $dval) = @$pair;
          $density_value += $dval * ($dlen/$total_len);
        }
      }
    }

    push @out, Bio::EnsEMBL::DensityFeature->new
      (-seq_region    => $slice,
       -start         => $start,
       -end           => $end,
       -density_type  => $density_value_type,
       -density_value => $density_value);

    $start = $end + 1;
    $end  += $wanted_block_size;
  }

  return \@out;
}


sub _tables {
  my $self = shift;

  return (['density_feature', 'df'], ['density_type', 'dt']);
}


sub _columns {
  my $self = shift;

  return qw( df.density_feature_id
             df.seq_region_id df.seq_region_start df.seq_region_end
             df.density_value dt.analysis_id dt.value_type);
}


sub _default_where_clause {
  my $self = shift;

  return 'df.density_type_id = dt.density_type_id';
}


sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches,
  # and a fair bit of gymnastics have been used.
  #

  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db()->get_AnalysisAdaptor();

  my @features;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my($density_feature_id, $seq_region_id, $seq_region_start, $seq_region_end,
     $density_value, $analysis_id, $density_value_type );

  $sth->bind_columns(\$density_feature_id, \$seq_region_id, \$seq_region_start,
                     \$seq_region_end, \$density_value, \$analysis_id, 
                     \$density_value_type);

  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;
  if($mapper) {
    $asm_cs = $mapper->assembled_CoordSystem();
    $cmp_cs = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $cmp_cs_vers = $cmp_cs->version();
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
  }

  FEATURE: while($sth->fetch()) {
    #get the analysis object
    my $analysis = $analysis_hash{$analysis_id} ||=
      $aa->fetch_by_dbID($analysis_id);

    #get the slice object
    my $slice = $slice_hash{"ID:".$seq_region_id};

    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    #
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    #
    if($mapper) {
      my $sr_name = $sr_name_hash{$seq_region_id};
      my $sr_cs   = $sr_cs_hash{$seq_region_id};

      my $len = $seq_region_end - $seq_region_start + 1;

      my @coords =
        $mapper->map($sr_name, $seq_region_start, $seq_region_end,1, $sr_cs);

      #filter out gaps
      @coords = grep {!$_->isa('Bio::EnsEMBL::Mapper::Gap')} @coords;

      #throw out density features mapped to gaps, or split
      next FEATURE if(@coords != 1);

      $seq_region_start  = $coords[0]->{'start'};
      $seq_region_end    = $coords[0]->{'end'};
      $sr_name           = $coords[0]->{'id'};

      if($density_value_type eq 'sum') {
        #adjust density value so it reflects length of feature actually used
        my $newlen = $seq_region_end - $seq_region_start +1;
        $density_value *= $newlen/$len if($newlen != $len);
      }

      #get a slice in the coord system we just mapped to
      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"} ||=
          $sa->fetch_by_region($cmp_cs_name, $sr_name,undef, undef, undef,
                               $cmp_cs_vers);
      } else {
        $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
          $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
                               $asm_cs_vers);
      }
    }

    #
    # If a destination slice was provided convert the coords
    # If the dest_slice starts at 1 and is foward strand, nothing needs doing
    #
    if($dest_slice) {
      if($dest_slice_start != 1 || $dest_slice_strand != 1) {
        if($dest_slice_strand == 1) {
          $seq_region_start = $seq_region_start - $dest_slice_start + 1;
          $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
        } else {
          my $tmp_seq_region_start = $seq_region_start;
          $seq_region_start = $dest_slice_end - $seq_region_end + 1;
          $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
        }

        #throw away features entirely off the end of the requested slice
        if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
          next FEATURE;
        }
      }
      $slice = $dest_slice;
    }

    push @features, Bio::EnsEMBL::DensityFeature->new
      (-start    => $seq_region_start,
       -end      => $seq_region_end,
       -seq_region    => $slice,
       #	-adaptor  => $self,
       -density_value      => $density_value,
       -density_type => $density_value_type); 
  }
  return \@features;
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$density_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all density features in the
               current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("density_feature");
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::DensityFeatures @df
               the simple features to store in the database
  Example    : $density_feature_adaptor->store(1234, @density_feats);
  Description: Stores a list of density feature objects in the database
  Returntype : none
  Exceptions : thrown if @df is not defined, if any of the features do not
               have an attached slice.
               or if any elements of @df are not Bio::EnsEMBL::SeqFeatures 
  Caller     : general

=cut

sub store{
  my ($self,@df) = @_;

  if( scalar(@df) == 0 ) {
    throw("Must call store with list of DensityFeatures");
  }
#mysql> desc density_feature;
#+--------------------+---------+------+-----+---------+----------------+
#| Field              | Type    | Null | Key | Default | Extra          |
#+--------------------+---------+------+-----+---------+----------------+
#| density_feature_id | int(11) |      | PRI | NULL    | auto_increment |
#| density_type_id    | int(11) |      | MUL | 0       |                |
#| seq_region_id      | int(11) |      |     | 0       |                |
#| seq_region_start   | int(11) |      |     | 0       |                |
#| seq_region_end     | int(11) |      |     | 0       |                |
#| density_value      | float   |      |     | 0       |                |
#+--------------------+---------+------+-----+---------+----------------+

  my $sth = $self->prepare
    ("INSERT INTO density_feature (seq_region_id, seq_region_start, " .
                                  "seq_region_end, density_type_id, " .
                                  "density_value) " .
     "VALUES (?,?,?,?,?)");

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

 FEATURE: foreach my $df ( @df ) {

    if( !ref $df || !$df->isa("Bio::EnsEMBL::DensityFeature") ) {
      throw("DensityFeature must be an Ensembl DensityFeature, " .
            "not a [".ref($df)."]");
    }

    if($df->is_stored($db)) {
      warning("DensityFeature [".$df->dbID."] is already stored" .
              " in this database.");
      next FEATURE;
    }

    if(!defined($df->density_type)) {
      throw("A density type must be attached to the features to be stored.");
    }

    #store the density_type if it has not been stored yet
    
    if(!$df->density_type->is_stored($db)) {
      $df->density_type->store($df->density_type_id);
    }

    my $original = $df;
    my $seq_region_id;
    ($df, $seq_region_id) = $self->_pre_store($df);

    $sth->execute($seq_region_id, $df->start, $df->end,
                  $df->density_type->dbID, $df->density_value);

    $original->dbID($sth->{'mysql_insertid'});
    $original->adaptor($self);
  }
}

sub fetch_all_by_slice_and_density_type{
  my ($self,$slice,$density_type) = @_;

  my $sth = $self->prepare
    ("SELECT density_feature_id, seq_region_start, seq_region_end, density_value  ".
     "FROM density_feature ".
     "WHERE seq_region_id=".$slice->get_seq_region_id." and ".
     "density_type_id = ".$density_type->dbID);
  
  my @out=();
  my ($id, $start, $end, $value);
  $sth->execute();
  $sth->bind_columns(\$id, \$start, \$end, \$value);


  while($sth->fetch()) {

    push @out, Bio::EnsEMBL::DensityFeature->new
      (-start    => $start,
       -end      => $end,
       -seq_region    => $slice,
       #	-adaptor  => $self,
       -density_value      => $value,
       -density_type => $density_type); 
}
  
  return \@out;
}

 


1;

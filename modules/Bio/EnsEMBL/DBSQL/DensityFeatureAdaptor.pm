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

my $density_feature_adaptor = $database_adaptor->get_DensityFeatureAdaptor();
@density_features = @{$density_feature_adaptor->fetch_all_by_Slice($slice)};

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

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


sub fetch_all_by_Slice {
  my ($self, $slice, $logic_name, $num_blocks, $interpolate, $max_ratio) = @_;

  $num_blocks ||= 50;
  my $length = $slice->length();
  my $wanted_block_size = int($length/$num_blocks);

  my $analysis_adaptor = $self->db()->get_AnalysisAdaptor();
  my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

  if(!$analysis) {
    warning("No Analysis with logic_name=[$logic_name] found. " .
            'Returning empty list.');
    return [];
  }

  my $sth = $self->prepare("SELECT density_type_id, block_size, analysis_id  ".
                           "FROM   density_type " .
                           "WHERE  analysis_id = ?");
  $sth->execute($analysis->dbID());

  my $best_ratio = undef;
  my $density_type_id;
  my $row;

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
    }
  }

  #the ratio was not good enough, or this logic name was not in the
  #density type table, return an empty list
  if(!defined($best_ratio) ||
     (defined($max_ratio) && $best_ratio > $max_ratio)) {
    return [];
  }

  my $constraint = "dt.density_type_id = $density_type_id";

  my @features =
    @{$self->fetch_all_by_Slice_constraint($slice,$constraint)};

  #we don't want to interpolate if the ratio was very close
  $interpolate = 0 if($best_ratio < 1.05);

  return \@features if(!$interpolate);

  #interpolate the features into new features of a different size
  my @out;
  #sort the features on start position
  @features = sort({$a->start() <=> $b->end()} @features);

  #resize the features that were returned
  my $start = 1;
  my $end   = $start+$wanted_block_size;

  while($start < $length) {
    $end = $length if($end > $length);

    my $density_value = 0.0;
    my ($f, $fstart, $fend, $portion);

    #construct a new feature using all of the old density features that
    #we overlapped
    while($f = shift(@features) && $end > $f->{'start'}) {

      #what portion of this feature are we using to construct our new block?
      $fstart = ($f->{'start'} < $start) ? $start : $f->{'start'};
      $fend   = ($f->{'end'}   > $end  ) ? $end   : $f->{'end'};
      $fend   = $length if($fend > $length);
      $portion = ($fend-$fstart+1)/$f->length();

      #take a percentage of density value, depending on how much of the
      #feature we overlapped
      $density_value += $portion * $f->{'density_value'};

      #do not want to look at next feature if we only used part of this one:
      last if($fend < $f->{'end'});
    }

    #if we did not completely overlap the last feature, put it back on so
    #it can be partially used by the next block
    if($fend < $f->{'end'}) {
      unshift(@features, $f);
    }

    push @out, Bio::EnsEMBL::DensityFeature->new_fast
      ({'start'    => $start,
        'end'      => $end,
        'slice'    => $slice,
        'analysis' => $analysis,
        'adaptor'  => $self,
        'density_value' => $density_value});

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
             df.density_value dt.analysis_id );
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
     $density_value, $analysis_id );

  $sth->bind_columns(\$density_feature_id, \$seq_region_id, \$seq_region_start,
                     \$seq_region_end, \$density_value, \$analysis_id);

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

      ($sr_name,$seq_region_start,$seq_region_end) =
        $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
                         1, $sr_cs);

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($sr_name));

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

        #throw away features off the end of the requested slice
        if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
          next FEATURE;
        }
      }
      $slice = $dest_slice;
    }

    push @features, Bio::EnsEMBL::DensityFeature->new_fast(
      {'start'    => $seq_region_start,
       'end'      => $seq_region_end,
       'slice'    => $slice,
       'analysis' => $analysis,
       'adaptor'  => $self,
       'dbID'     => $density_feature_id,
       'density_value' => $density_value,});
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

1;

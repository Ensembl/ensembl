#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor

=head1 SYNOPSIS

$mfa = $database_adaptor->get_MiscFeatureAdaptor();

my $misc_feat = $mfa->fetch_by_dbID(1234);
my @misc_feats = @{$mfa->fetch_all_by_Slice($slice)};

my @misc_clones = @{$mfa->fetch_all_by_Slice_and_set_code('cloneset')};


=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of MiscFeatures.  Misc
Features are extremely generic features that can be added with minimal
effort to the database.  Currently misc features are used to describe
the locations of clone sets and tiling path information, but arbitrary
features can be stored.  Misc features are grouped into sets and can be
fetched according to their grouping using the fetch_all_by_Slice_and_set_code
and fetch_all_by_set_code methods.  MiscFeatures may belong to more
than one set.


=head1 CONTACT

Post questions to the EnsEMBL developer list: ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);



=head2 fetch_all_by_Slice_and_set_code

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               A slice representing the region to fetch from
  Arg [2]    : string $set_code
               The code of the set to retrieve features from
  Example    : @feats = @{$mfa->fetch_all_by_Slice_and_set_code('cloneset')};
  Description: Retrieves a set of MiscFeatures which have a particular set code
               and which lie in a particular region.  All features with the
               provide set code and which overlap the given slice are returned.
  Returntype : listref of Bio::EnsEMBL::MiscFeatures
  Exceptions : throw if set_code is not provided
               warning if no set for provided set code exists
  Caller     : general

=cut

sub fetch_all_by_Slice_and_set_code {
  my $self = shift;
  my $slice = shift;
  my $set_code = shift;

  throw('Set code argument is required.') if(!$set_code);

  my $msa = $self->db->get_MiscSetAdaptor();
  my $set = $msa->fetch_by_code($set_code);

  if(!$set) {
    warning("No misc_set with code [$set_code] exists.\n" .
            "Returning empty list.");
    return [];
  }

  my $constraint = " mfms.misc_set_id = " . $set->dbID();
  return $self->fetch_all_by_Slice_constraint($slice, $constraint);
}



=head2 fetch_all_by_attribute_type_value

  Arg [1]    : string $attrib_type_code
               The code of the attribute type to fetch features for
  Arg [2]    : (optional) string $attrib_value
               The value of the attribute to fetch features for
  Example    : 
         #get all misc features that have an embl accession
         @feats = @{$mfa->fetch_all_by_attrib_type_value('embl_acc')};
         #get the misc feature with synonym 'AL014121'
         ($feat)=@{$mfa->fetch_all_by_attrib_type_value('synonym','AL014121');
  Description: Retrieves MiscFeatures which have a particular attribute.
               If the attribute value argument is also provided only features
               which have the attribute AND a particular value are returned.
               The features are returned in their native coordinate system
               (i.e. the coordinate system that they are stored in).
  Returntype : listref of Bio::EnsEMBL::MiscFeatures
  Exceptions : throw if attrib_type code arg is not provided
  Caller     : general

=cut

sub fetch_all_by_attribute_type_value {
  my $self = shift;
  my $attrib_type_code = shift;
  my $attrib_value = shift;

  throw("Attrib type code argument is required.") if(!$attrib_type_code);

  #
  # Need to do 2 queries so that all of the ids come back with the features.
  # The problem with adding attrib constraints to filter the misc_features
  # which come back is that not all of the attributes will come back
  #

  my $sql = "SELECT DISTINCT ma.misc_feature_id " .
            "FROM   misc_attrib ma, attrib_type at " .
            "WHERE  ma.attrib_type_id = at.attrib_type_id " .
            "AND    at.code = ?";

  my @bind_vals = ($attrib_type_code);

  if($attrib_value) {
    push @bind_vals, $attrib_value;
    $sql .= " AND ma.value = ?";
  }

  my $sth = $self->prepare($sql);
  $sth->execute(@bind_vals);

  my @ids = map {$_->[0]} @{$sth->fetchall_arrayref()};

  $sth->finish();

  #construct constraints from the list of ids.  Split ids into
  #groups of 1000 to ensure that the query is not too big
  my @constraints;
  while(@ids) {
    my @subset =  splice(@ids, 0, 1000);
    if(@subset == 1) {
      push @constraints, "mf.misc_feature_id = $subset[0]";
    } else {
      my $id_str = join(',',@subset);
      push @constraints, "mf.misc_feature_id in ($id_str)";
    }
  }

  my @results;
  foreach my $constraint (@constraints) {
    push @results, @{$self->generic_fetch($constraint)};
  }

  return \@results;
}


#_tables
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED Implementation of abstract superclass method to
#               provide the name of the tables to query
#  Returntype : string
#  Exceptions : none
#  Caller     : internal


sub _tables {
  my $self = shift;

  return (['misc_feature',          'mf'],
          ['misc_feature_misc_set', 'mfms'],
          ['misc_attrib',           'ma'],
          ['attrib_type',           'at']);
}


#_columns

#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED Implementation of abstract superclass method to 
#               provide the name of the columns to query 
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal

sub _columns {
  my $self = shift;

  #warning _objs_from_sth implementation depends on ordering
  return qw (mf.misc_feature_id
	     mf.seq_region_id
	     mf.seq_region_start
	     mf.seq_region_end
	     mf.seq_region_strand
	     ma.value
       at.code
       mfms.misc_set_id);
}



# _default_where_clause

#  Arg [1]    : none
#  Example    : none
#  Description: Overrides superclass method to provide an additional 
#               table joining constraint before the SQL query is performed.
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch

sub _default_where_clause {
  my $self = shift;

  return '';
}


sub _left_join {
  my $self = shift;

  return(
      ['misc_feature_misc_set','mf.misc_feature_id = mfms.misc_feature_id'],
      ['misc_attrib', 'mf.misc_feature_id = ma.misc_feature_id'],
      ['attrib_type','ma.attrib_type_id = at.attrib_type_id']);
}


sub _final_clause {
  my $self = shift;

  return "order by mf.misc_feature_id";
}


# _objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Example    : none
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of MiscFeatures from a
#               hashref generated from an SQL query
#  Returntype : listref of Bio::EnsEMBL::MiscFeatures
#  Exceptions : none
#  Caller     : internal

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $sa = $self->db()->get_SliceAdaptor();
  my $msa = $self->db->get_MiscSetAdaptor();

  my @features;
  my %ms_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my($misc_feature_id, $seq_region_id, $seq_region_start, $seq_region_end,
     $seq_region_strand, $attrib_value, $attrib_type_code, $misc_set_id);

  $sth->bind_columns( \$misc_feature_id, \$seq_region_id, \$seq_region_start,
                      \$seq_region_end, \$seq_region_strand,
                      \$attrib_value, \$attrib_type_code,\$misc_set_id);

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

  my $current = -1;
  my $throw_away = -1;
  my $feat;
  my $feat_misc_sets;
  my $feat_attribs;
  my $seen_attribs;

 FEATURE: while($sth->fetch()) {
    #if this feature is not being used, skip all rows related to it
    next if($throw_away == $misc_feature_id);

    if($current == $misc_feature_id) {
      #still working on building up attributes and sets for current feature

      #if there is a misc_set, add it to the current feature
      if($misc_set_id) {
        my $misc_set = $ms_hash{$misc_set_id} ||=
          $msa->fetch_by_dbID($misc_set_id);
        #doesn't matter if it is added twice in same slot
        $feat_misc_sets->{$misc_set->{'code'}} = $misc_set;
      }

      #if there is a new attribute add it to the current feature
      if($attrib_value && $attrib_type_code &&
         !$seen_attribs->{"$attrib_type_code:$attrib_value"}) {
        $feat_attribs->{$attrib_type_code} ||= [];
        push @{$feat_attribs->{$attrib_type_code}}, $attrib_value;
        $seen_attribs->{"$attrib_type_code:$attrib_value"} = 1;
      }

    } else {
      if($feat) {
        #start working on a new feature, discard references to last one
        $feat = {};
        $feat_attribs = {};
        $feat_misc_sets = {};
        $seen_attribs = {};
      }

      $current = $misc_feature_id;
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

        ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
          $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
                           $seq_region_strand, $sr_cs);

        #skip features that map to gaps or coord system boundaries
        next FEATURE if(!defined($sr_name));

        #get a slice in the coord system we just mapped to
        if($asm_cs == $sr_cs ||
           ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
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
      if($dest_slice && ($dest_slice_start != 1 || $dest_slice_strand != 1)) {
        if($dest_slice_strand == 1) {
          $seq_region_start = $seq_region_start - $dest_slice_start + 1;
          $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
        } else {
          my $tmp_seq_region_start = $seq_region_start;
          $seq_region_start = $dest_slice_end - $seq_region_end + 1;
          $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
          $seq_region_strand *= -1;
        }

        $slice = $dest_slice;

        #throw away features off the end of the requested slice
        if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
          #flag this feature as one to throw away
          $throw_away = $misc_feature_id;
          next FEATURE;
        }
      }

      if($misc_set_id) {
        #get the misc_set object
        my $misc_set = $ms_hash{$misc_set_id} ||=
          $msa->fetch_by_dbID($misc_set_id);
        $feat_misc_sets->{$misc_set->{'code'}} = $misc_set;
      }

      if($attrib_value && $attrib_type_code) {
        $feat_attribs->{$attrib_type_code} = [$attrib_value];
        $seen_attribs->{"$attrib_type_code:$attrib_value"} = 1;
      }

      $feat =  Bio::EnsEMBL::MiscFeature->new_fast
        ({'start'   => $seq_region_start,
          'end'     => $seq_region_end,
          'strand'  => $seq_region_strand,
          'slice'   => $slice,
          'adaptor' => $self,
          'dbID'    => $misc_feature_id,
          'attributes' => $feat_attribs,
          'sets'    => $feat_misc_sets});
      push @features, $feat;
    }
  }

  return \@features;
}



=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$misc_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all misc_features in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("misc_feature");
}

1;






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

Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor

=head1 SYNOPSIS

  $mfa = $database_adaptor->get_MiscFeatureAdaptor();

  # retrieve a misc feature by its dbID
  my $misc_feat = $mfa->fetch_by_dbID(1234);

  # retrieve all misc features in a given region
  my @misc_feats = @{ $mfa->fetch_all_by_Slice($slice) };

  # retrieve all misc features in a given region with a given set code
  my @misc_clones =
    @{ $mfa->fetch_all_by_Slice_and_set_code('cloneset') };

  # store some misc features in the database
  $mfa->store(@misc_features);

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of MiscFeatures.
Misc Features are extremely generic features that can be added with
minimal effort to the database.  Currently misc features are used to
describe the locations of clone sets and tiling path information,
but arbitrary features can be stored.  Misc features are grouped
into sets and can be fetched according to their grouping using the
fetch_all_by_Slice_and_set_code and fetch_all_by_set_code methods.
MiscFeatures may belong to more than one set.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);



=head2 fetch_all_by_Slice_and_set_code

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               A slice representing the region to fetch from
  Arg [2...] : string $set_code
               The code of the set to retrieve features from
  Example    : @feats = @{$mfa->fetch_all_by_Slice_and_set_code('cloneset')};
  Description: Retrieves a set of MiscFeatures which have a particular set code
               and which lie in a particular region.  All features with the
               provide set code and which overlap the given slice are returned.
  Returntype : listref of Bio::EnsEMBL::MiscFeatures
  Exceptions : throw if set_code is not provided
               warning if no set for provided set code exists
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_and_set_code {
  my $self = shift;
  my $slice = shift;

  throw('Set code argument is required.') unless @_;

  my $msa = $self->db->get_MiscSetAdaptor();
  my @sets = ();
  my $max_len = 0;
  foreach my $set_code (@_) {
    my $set = $msa->fetch_by_code($set_code);
    if($set) {
      $max_len = $set->longest_feature if $set->longest_feature > $max_len;
      push @sets, $set->dbID;
    } else { 
      warning("No misc_set with code [$set_code] exists");
    }
  }
  my $constraint;
  if( @sets > 1 ) {
    $constraint = " mfms.misc_set_id in ( @{[join ',',@sets]} ) ";
  } elsif( @sets == 1 ) {
    $constraint = " mfms.misc_set_id = $sets[0] ";
  } else {
    return [];
  }

  $self->_max_feature_length($max_len);

  my $results = $self->fetch_all_by_Slice_constraint($slice, $constraint);

  $self->_max_feature_length(undef);

  return $results;
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
               If the attribute value argument is also provided only
               features which have the attribute AND a particular value
               are returned.  The features are returned in their native
               coordinate system (i.e. the coordinate system that they
               are stored in).
  Returntype : listref of Bio::EnsEMBL::MiscFeatures
  Exceptions : throw if attrib_type code arg is not provided
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_attribute_type_value {
  my $self             = shift;
  my $attrib_type_code = shift;
  my $attrib_value     = shift;

  throw("Attrib type code argument is required.")
    if ( !$attrib_type_code );

  # Need to do 2 queries so that all of the ids come back with the
  # features.  The problem with adding attrib constraints to filter the
  # misc_features which come back is that not all of the attributes will
  # come back

  my $sql = qq(
  SELECT DISTINCT
        ma.misc_feature_id
  FROM  misc_attrib ma,
        attrib_type at,
        misc_feature mf,
        seq_region sr,
        coord_system cs
  WHERE ma.attrib_type_id = at.attrib_type_id
    AND at.code = ?
    AND ma.misc_feature_id = mf.misc_feature_id
    AND mf.seq_region_id = sr.seq_region_id
    AND sr.coord_system_id = cs.coord_system_id
    AND cs.species_id = ?);

  if ($attrib_value) {
    $sql .= " AND ma.value = ?";
  }

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $attrib_type_code,   SQL_VARCHAR );
  $sth->bind_param( 2, $self->species_id(), SQL_INTEGER );
  if ($attrib_value) {
    $sth->bind_param( 3, $attrib_value, SQL_VARCHAR );
  }

  $sth->execute();

  my @ids = map { $_->[0] } @{ $sth->fetchall_arrayref() };

  $sth->finish();

  # Construct constraints from the list of ids.  Split ids into groups
  # of 1000 to ensure that the query is not too big.
  my @constraints;
  while (@ids) {
    my @subset = splice( @ids, 0, 1000 );
    if ( @subset == 1 ) {
      push @constraints, "mf.misc_feature_id = $subset[0]";
    } else {
      my $id_str = join( ',', @subset );
      push @constraints, "mf.misc_feature_id in ($id_str)";
    }
  }

  my @results;
  foreach my $constraint (@constraints) {
    push @results, @{ $self->generic_fetch($constraint) };
  }

  return \@results;
} ## end sub fetch_all_by_attribute_type_value


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
	     mfms.misc_set_id
	     at.name
	     at.description);
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

  return " ORDER BY mf.misc_feature_id";
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
     $seq_region_strand, $attrib_value, $attrib_type_code, $misc_set_id,
     $attrib_type_name, $attrib_type_description );

  $sth->bind_columns( \$misc_feature_id, \$seq_region_id, \$seq_region_start,
                      \$seq_region_end, \$seq_region_strand,
                      \$attrib_value, \$attrib_type_code,\$misc_set_id,
		      \$attrib_type_name, \$attrib_type_description );

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
  my $dest_slice_sr_name;
  my $dest_slice_sr_id;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
    $dest_slice_sr_id  = $dest_slice->get_seq_region_id();
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

    if ($current == $misc_feature_id) {
      #still working on building up attributes and sets for current feature

      #if there is a misc_set, add it to the current feature
      if ($misc_set_id) {
        my $misc_set = $ms_hash{$misc_set_id} ||=
          $msa->fetch_by_dbID($misc_set_id);
        if ( ! exists $feat_misc_sets->{$misc_set->{'code'}} ) {
          $feat->add_MiscSet( $misc_set );
          $feat_misc_sets->{$misc_set->{'code'}} = $misc_set;
        }
      }

      #if there is a new attribute add it to the current feature
      if ($attrib_value && $attrib_type_code &&
          !$seen_attribs->{"$attrib_type_code:$attrib_value"}) {
        my $attrib = Bio::EnsEMBL::Attribute->new
          ( -CODE => $attrib_type_code,
            -NAME => $attrib_type_name,
            -DESC => $attrib_type_description,
            -VALUE => $attrib_value
          );
	
	
        $feat_attribs ||= [];
        push @$feat_attribs, $attrib;
        $seen_attribs->{"$attrib_type_code:$attrib_value"} = 1;
      }

    } else {
      if ($feat) {
        #start working on a new feature, discard references to last one
        $feat = {};
        $feat_attribs = [];
        $feat_misc_sets = {};
        $seen_attribs = {};
      }

      $current = $misc_feature_id;
      #need to get the internal_seq_region, if present
      $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
      my $slice = $slice_hash{"ID:".$seq_region_id};

      if (!$slice) {
        $slice = $sa->fetch_by_seq_region_id($seq_region_id);
        $slice_hash{"ID:".$seq_region_id} = $slice;
        $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
        $sr_cs_hash{$seq_region_id} = $slice->coord_system();
      }

      my $sr_name = $sr_name_hash{$seq_region_id};
      my $sr_cs   = $sr_cs_hash{$seq_region_id};
      #
      # remap the feature coordinates to another coord system
      # if a mapper was provided
      #
      if ($mapper) {


	if (defined $dest_slice && $mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper')  ) {
	    ( $seq_region_id,  $seq_region_start,$seq_region_end, $seq_region_strand ) =
	       $mapper->map( $sr_name, $seq_region_start, $seq_region_end,
                          $seq_region_strand, $sr_cs, 1, $dest_slice);

	} else {

	    ( $seq_region_id,  $seq_region_start,$seq_region_end, $seq_region_strand ) =
	       $mapper->fastmap( $sr_name, $seq_region_start, $seq_region_end,$seq_region_strand, $sr_cs );
	}

        #skip features that map to gaps or coord system boundaries
        if(!defined($seq_region_id)) {
          $throw_away = $misc_feature_id;
          next FEATURE;
        }

        #get a slice in the coord system we just mapped to
#        if ($asm_cs == $sr_cs ||
#            ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
          $slice = $slice_hash{"ID:".$seq_region_id} ||=
            $sa->fetch_by_seq_region_id($seq_region_id);
#        } else {
#          $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
#            $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
#                                 $asm_cs_vers);
#        }
      }

      #
      # If a destination slice was provided convert the coords
      # If the dest_slice starts at 1 and is foward strand, nothing needs doing
      #
      if ($dest_slice) {
        if ($dest_slice_start != 1 || $dest_slice_strand != 1) {
          if ($dest_slice_strand == 1) {
            $seq_region_start = $seq_region_start - $dest_slice_start + 1;
            $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
          } else {
            my $tmp_seq_region_start = $seq_region_start;
            $seq_region_start = $dest_slice_end - $seq_region_end + 1;
            $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
            $seq_region_strand *= -1;
          }
	}
	#throw away features off the end of the requested slice
	if ($seq_region_end < 1 || $seq_region_start > $dest_slice_length ||
	   ( $dest_slice_sr_id ne $seq_region_id )) {
	  #flag this feature as one to throw away
	  $throw_away = $misc_feature_id;
	  next FEATURE;
	}
        $slice = $dest_slice;
      }


      if ($attrib_value && $attrib_type_code) {
        my $attrib = Bio::EnsEMBL::Attribute->new
          ( -CODE => $attrib_type_code,
            -NAME => $attrib_type_name,
            -DESC => $attrib_type_description,
            -VALUE => $attrib_value
          );
        $feat_attribs = [$attrib];
        $seen_attribs->{"$attrib_type_code:$attrib_value"} = 1;
      }

      $feat =
        $self->_create_feature_fast( 'Bio::EnsEMBL::MiscFeature', {
                                       'start'   => $seq_region_start,
                                       'end'     => $seq_region_end,
                                       'strand'  => $seq_region_strand,
                                       'slice'   => $slice,
                                       'adaptor' => $self,
                                       'dbID'    => $misc_feature_id,
                                       'attributes' => $feat_attribs
                                     } );

      push @features, $feat;

      if ($misc_set_id) {
        #get the misc_set object
        my $misc_set = $ms_hash{$misc_set_id} ||=
          $msa->fetch_by_dbID($misc_set_id);
        if ( ! exists $feat_misc_sets->{$misc_set->{'code'}} ) {
          $feat->add_MiscSet( $misc_set );
          $feat_misc_sets->{$misc_set->{'code'}} = $misc_set;
        }
      }
    }
  }

  return \@features;
}



=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$misc_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all misc_features in the 
               current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self,$ordered) = @_;

   return $self->_list_dbIDs("misc_feature",undef,$ordered);
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::MiscFeatures @misc_features
  Example    : $misc_feature_adaptor->store(@misc_features);
  Description: Stores a list of MiscFeatures in this database.  The stored
               features will have their 
  Returntype : none
  Exceptions : throw on invalid arguments
               warning if misc feature is already stored in this database
               throw if start/end/strand attribs are not valid
  Caller     : general
  Status     : Stable

=cut

sub store {
  my $self = shift;
  my @misc_features = @_;

  my $db = $self->db();

  my $feature_sth = $self->prepare
    ("INSERT INTO misc_feature SET " .
     " seq_region_id    = ?, " .
     " seq_region_start = ?, " .
     " seq_region_end   = ?, " .
     " seq_region_strand = ?");

  my $feature_set_sth = $self->prepare
    ("INSERT IGNORE misc_feature_misc_set SET " .
     " misc_feature_id = ?, " .
     " misc_set_id = ?");

  my $msa = $db->get_MiscSetAdaptor();
  my $aa  = $db->get_AttributeAdaptor();

 FEATURE:
  foreach my $mf (@misc_features) {
    if(!ref($mf) || !$mf->isa('Bio::EnsEMBL::MiscFeature')) {
      throw("List of MiscFeature arguments expeceted");
    }

    if($mf->is_stored($db)) {
      warning("MiscFeature [" .$mf->dbID."] is already stored in database.");
      next FEATURE;
    }

    # do some checking of the start/end and convert to seq_region coords
    my $original = $mf;
    my $seq_region_id;
    ($mf, $seq_region_id) = $self->_pre_store($mf);

    # store the actual MiscFeature
    $feature_sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $feature_sth->bind_param(2,$mf->start,SQL_INTEGER);
    $feature_sth->bind_param(3,$mf->end,SQL_INTEGER);
    $feature_sth->bind_param(4,$mf->strand,SQL_TINYINT);
    $feature_sth->execute();

    my $dbID = $feature_sth->{'mysql_insertid'};

    $mf->dbID($dbID);
    $mf->adaptor($self);

    # store all the attributes
    my $attribs = $mf->get_all_Attributes();
    $aa->store_on_MiscFeature($mf, $attribs);

    # store all the sets that have not been stored yet
    my $sets = $mf->get_all_MiscSets();
    foreach my $set (@$sets) {
      $msa->store($set) if(!$set->is_stored($db));

      # update the misc_feat_misc_set table to store the set relationship
      $feature_set_sth->bind_param(1,$dbID,SQL_INTEGER);
      $feature_set_sth->bind_param(2,$set->dbID,SQL_INTEGER);

      $feature_set_sth->execute();
    }
  }

  return;
}

1;






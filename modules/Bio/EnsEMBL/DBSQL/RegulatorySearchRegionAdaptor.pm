#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::RegulatorySearchRegionAdaptor
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::RegulatorySearchRegionAdaptor

=head1 SYNOPSIS

$rfa = $database_adaptor->get_RegulatorySearchRegionAdaptor();

my $rs = $rfa->fetch_by_dbID(1234);

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of RegulatorySearchRegion objects
from the database.  Most of the implementation is in the superclass BaseFeatureAdaptor.

=head1 AUTHOR - Glenn Proctor

=head1 CONTACT

Post questions to the EnsEMBL developer list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::RegulatorySearchRegionAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::RegulatorySearchRegion;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);




=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : (optional) string $logic_name
               Limits RegulatorySearchRegions obtained to those having an Analysis with
               of the specified logic_name.  If no logic name is specified,
               regulatory features of all analysis types are retrieved.
  Example    : @rss = @{$rsa->fetch_all_by_Slice($slice, 'CisRed_search')};
  Description: Retrieves regulatory search regions overlapping the area designated by
               the provided slice argument.  Returned regions will be in
               in the same coordinate system as the provided slice and will
               have coordinates relative to the slice start.
  Returntype : reference to a list of Bio::EnsEMBL::RegulatorySearchRegions.
  Exceptions : throw on bad argument
  Caller     : Slice::get_all_RegulatorySearchRegions
  Status     : At Risk
             : under development

=cut

sub fetch_all_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $logic_name = shift;

  my $result = $self->fetch_all_by_Slice_constraint($slice,undef,$logic_name);

  return $result;
}


#  _tablename
#
#   Arg [1]    : none
#   Example    : none
#   Description: PROTECTED Implementation of abstract superclass method to 
#                provide the name of the tables to query 
#   Returntype : string
#   Exceptions : none
#   Caller     : internal


sub _tables {
  my $self = shift;

  return (['regulatory_search_region', 'rs']);
}


# _columns
#
#   Arg [1]    : none
#   Example    : none
#   Description: PROTECTED Implementation of abstract superclass method to 
#                provide the name of the columns to query 
#   Returntype : list of strings
#   Exceptions : none
#   Caller     : internal

sub _columns {
  my $self = shift;

  return qw (rs.regulatory_search_region_id
	     rs.name
	     rs.seq_region_id
	     rs.seq_region_start
	     rs.seq_region_end
	     rs.seq_region_strand
	     rs.analysis_id
	     rs.ensembl_object_id
	     rs.ensembl_object_type);
}


# _default_where_clause
#  Arg [1]    : none
#  Example    : none
#  Description: Overrides superclass method to provide an additional
#               table joining constraint before the SQL query is performed.
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch
#

sub _default_where_clause {
  my $self = shift;

  return '';
}



#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of RegulatorySearchRegions from a
#               hashref generated from an SQL query

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db->get_AnalysisAdaptor();

  my @features;
  my %rm_hash;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my($regulatory_search_region_id, $name, $seq_region_id,
     $seq_region_start, $seq_region_end, $seq_region_strand, $analysis_id,
     $ensembl_object_id, $ensembl_object_type);

  $sth->bind_columns( \$regulatory_search_region_id, \$name, \$seq_region_id, 
		      \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$analysis_id,
		      \$ensembl_object_id, \$ensembl_object_type);

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
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
  }

  FEATURE: while($sth->fetch()) {
    # create a regulatory search region object

    #get the analysis object
    my $analysis = $analysis_hash{$analysis_id} ||=
      $aa->fetch_by_dbID($analysis_id);

    my $slice = $slice_hash{"ID:".$seq_region_id};

    if(!$slice) {
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
    if($mapper) {

      ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
        $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
                          $seq_region_strand, $sr_cs);

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
          $seq_region_strand *= -1;
        }
      }

      #throw away features off the end of the requested slice
      if($seq_region_end < 1 || $seq_region_start > $dest_slice_length ||
	( $dest_slice_sr_name ne $sr_name )) {
	next FEATURE;
      }
      $slice = $dest_slice;
    }

    #finally, create the new regulatory search region
    push @features, Bio::EnsEMBL::RegulatorySearchRegion->new_fast
      ( { 'analysis'                 =>  $analysis,
	  'name'                     =>  $name,
          'start'                    =>  $seq_region_start,
          'end'                      =>  $seq_region_end,
          'strand'                   =>  $seq_region_strand,
          'ensembl_object_id'        =>  $ensembl_object_id,
          'ensembl_object_type'      =>  $ensembl_object_type,
          'adaptor'                  =>  $self,
          'slice'                    =>  $slice,
          'dbID'                     =>  $regulatory_search_region_id } );

  }

  return \@features;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::RegulatorySearchRegion
               the regulatory feature to store in the database
  Example    : $regulatory_search_region_adaptor->store($regulatory_search_region);
  Description: stores regulatory features in the database
  Returntype : none
  Exceptions :
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub store {
  my( $self, $feature, $ensObj, $ensType ) = @_;

  my $rs_sth = $self->prepare(qq {INSERT into regulatory_search_region
				  (name,
				   seq_region_id,
				   seq_region_start,
				   seq_region_end,
				   seq_region_strand,
				   analysis_id,
				   ensembl_object_type,
				   ensembl_object_id)
				  VALUES (?,?,?,?,?,?,?,?)});

  if (!ref($feature) || !$feature->isa('Bio::EnsEMBL::RegulatorySearchRegion')) {
    throw('Expected RegulatorySearchRegion argument not [' . ref($feature) .'].');
  }

  my $name = $feature->name() or throw("name not set");

  my $analysis = $feature->analysis();
  if (!ref($analysis) || !$analysis->isa("Bio::EnsEMBL::Analysis")) {
    throw("RegulatorySearchRegion cannot be stored without an associated analysis.");
  }

  my $original = $feature;
  my $seq_region_id;
  ($feature, $seq_region_id) = $self->_pre_store($feature);

  $rs_sth->execute($feature->name(),
		   $seq_region_id,
		   $feature->start(),
		   $feature->end(),
		   $feature->strand(),
		   $analysis->dbID(),
		   $feature->ensembl_object_type(),
		   $feature->ensembl_object_id());

  my $db_id = $rs_sth->{'mysql_insertid'}
    or throw("Didn't get an insertid from the INSERT statement");

  $original->dbID($db_id);
  $original->adaptor($self);

}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$repeat_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all repeat features in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : At Risk
             : under development

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("regulatory_search_region");
}


=head2 fetch_by_name

  Arg [1]    : string $name
               the name of the regulatory search_region to obtain
  Example    : $rs = $rsa->fetch_by_name('CisRed_Search_11');
  Description: Obtains a regulatory factor from the database via its name
  Returntype : Bio::EnsEMBL::RegulatorySearchRegion
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub fetch_by_name {
    my( $self, $name ) = @_;

    my ($rc) = @{$self->generic_fetch("name = '$name'")};

    return $rc;
}


1;



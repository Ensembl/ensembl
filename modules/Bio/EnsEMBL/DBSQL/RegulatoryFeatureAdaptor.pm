#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::RegulatoryFeatureAdaptor
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::RegulatoryFeatureAdaptor

=head1 SYNOPSIS

$rfa = $database_adaptor->get_RegulatoryFeatureAdaptor();

my $regulatory_feature = $rfa->fetch_by_dbID(1234);

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of RegulatoryFeature objects
from the database.  Most of the implementation is in the superclass BaseFeatureAdaptor.

=head1 AUTHOR - Glenn Proctor

=head1 CONTACT

Post questions to the EnsEMBL developer list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::RegulatoryFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::RegulatoryFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);




=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : (optional) string $logic_name
               Limits RepeatFeatures obtained to those having an Analysis with
               of the specified logic_name.  If no logic name is specified
               Repeats of all analysis types are retrieved.
  Example    : @rfeats = @{$rfa->fetch_all_by_Slice($slice, 'miRanda')};
  Description: Retrieves regulatory features overlapping the area designated by
               the provided slice argument.  Returned features will be in
               in the same coordinate system as the provided slice and will
               have coordinates relative to the slice start.
  Returntype : reference to a list of Bio::EnsEMBL::RegulatoryFeatures.
  Exceptions : throw on bad argument
  Caller     : Slice::get_all_RegulatoryFeatures

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

  return (['regulatory_feature', 'rf'], ['regulatory_motif', 'rm']);
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

  return qw (rf.regulatory_feature_id
	     rf.name
	     rf.seq_region_id
	     rf.seq_region_start
	     rf.seq_region_end
	     rf.seq_region_strand
	     rf.analysis_id
	     rf.influence
	     rf.regulatory_motif_id
	     rm.name
	     rm.type);
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

  return 'rf.regulatory_motif_id = rm.regulatory_motif_id';
}



#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of RegulatoryFeatures from a
#               hashref generated from an SQL query

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $rca = $self->db()->get_RegulatoryMotifAdaptor();
  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db->get_AnalysisAdaptor();

  my @features;
  my %rm_hash;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my($regulatory_feature_id, $regulatory_feature_name, $seq_region_id,
     $seq_region_start, $seq_region_end, $seq_region_strand, $analysis_id,
     $influence, $motif_id, $motif_name, $motif_type);

  $sth->bind_columns( \$regulatory_feature_id, \$regulatory_feature_name, \$seq_region_id, 
		      \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$analysis_id,
		      \$influence, \$motif_id, \$motif_name, \$motif_type);

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
    #create a regulatory motif object

    my $rm = $rm_hash{$motif_id} ||=
      Bio::EnsEMBL::RegulatoryMotif->new_fast
          ({'dbID' => $motif_id,
            'name' => $motif_name,
            'type' => $motif_type});

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

    #finally, create the new regulatory feature
    push @features, Bio::EnsEMBL::RegulatoryFeature->new_fast
      ( { 'analysis'      =>  $analysis,
	  'name'          =>  $regulatory_feature_name,
          'start'         =>  $seq_region_start,
          'end'           =>  $seq_region_end,
          'strand'        =>  $seq_region_strand,
          'motif'         =>  $rm,
          'influence'     =>  $influence,
          'adaptor'       =>  $self,
          'slice'         =>  $slice,
          'dbID'          =>  $regulatory_feature_id } );

  }

  return \@features;
}

=head2 fetch_all_by_motif

  Arg [1]    : Bio::EnsEMBL::RegulatoryMotif 
               the type of regulatory motif to obtain
  Example    : $rm = $rma->fetch_all_by_motif($motif);
  Description: Obtains all regulatory features that correspond to a
               particular regulatory motif
  Returntype : listREF of Bio::EnsEMBL::RegulatoryFeatures
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_motif {
    my( $self, $motif) = @_;

    return $self->generic_fetch("rf.regulatory_motif_id = " . $motif->dbID());
}



=head2 store

  Arg [1]    : list of Bio::EnsEMBL::RegulatoryFeatures
               the regulatory features to store in the database
  Example    : $regulatory_feature_adaptor->store($regulatory_feature);
  Description: stores regulatory features in the database
  Returntype : none
  Exceptions :
  Caller     : general

=cut

sub store {
  my( $self, @features ) = @_;

  my $sth = $self->prepare(qq {INSERT into regulatory_feature 
                           (name,
                            seq_region_id,
                            seq_region_start,
                            seq_region_end,
                            seq_region_strand,
                            analysis_id,
                            regulatory_motif_id,
                            influence) 
                            VALUES (?,?,?,?,?,?,?,?)});

  foreach my $rf (@features) {

    if(!ref($rf) || !$rf->isa('Bio::EnsEMBL::RegulatoryFeature')) {
      throw('Expected RegulatoryFeature argument not [' . ref($rf) .'].');
    }

    my $name = $rf->name() or throw("name not set");

    my $analysis = $rf->analysis();
    if(!ref($analysis) || !$analysis->isa("Bio::EnsEMBL::Analysis")) {
      throw("RegulatoryFeature cannot be stored without an associated analysus.");
    }

    my $original = $rf;
    my $seq_region_id;
    ($rf, $seq_region_id) = $self->_pre_store($rf);

    $sth->execute($rf->name(),
		  $seq_region_id,
		  $rf->start(),
		  $rf->end(),
		  $rf->strand(),
		  $analysis->dbID(),
		  $rf->motif()->dbID(),
		  $rf->influence());

    my $db_id = $sth->{'mysql_insertid'}
    or throw("Didn't get an insertid from the INSERT statement");

    $original->dbID($db_id);
    $original->adaptor($self);
  }
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$repeat_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all repeat features in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("regulatory_feature");
}

=head2 fetch_all_by_ensembl_object_type

  Arg [1]    : string $type - one of 'Gene', 'Transcript', 'Translation'
  Arg [2]    : dbID of gene/transcript/translation
  Example    : @features = @{$regulatory_feature_adaptor->
                      fetch_all_by_ensembl_object_type('Transcript', 21050)};
  Description: Gets all the regulatory features associated with a particular
               gene, transcript or translation. Each feature only appears once.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature
  Exceptions : none
  Caller     : ?

=cut

sub fetch_all_by_ensembl_object_type {

   my ($self, $type, $id) = @_;

   my $sth = $self->prepare("SELECT regulatory_feature_id FROM regulatory_feature_object WHERE ensembl_object_type=? AND ensembl_object_id=?");

   $sth->execute($type, $id);

   my $dbID;
   my %regulatory_features;
   while (($dbID) = $sth->fetchrow_array()) {
     my $feature = $self->fetch_by_dbID($dbID);
     if (!exists($regulatory_features{$feature->dbID()})) {
       $regulatory_features{$feature->dbID()} = $feature;
     }
   }

   my @features = values %regulatory_features;

   return \@features;

}

=head2 fetch_all_by_transcript

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : @features = @{$regulatory_feature_adaptor->
                      fetch_all_by_transcript($transcript)};
  Description: Gets all the regulatory features associated with a
               particular transcript. Each feature only appears once.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature
  Exceptions : If arg is not of correct type.
  Caller     : ?

=cut

sub fetch_all_by_transcript {

   my ($self, $transcript) = @_;

   if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
     throw('Expected Bio::EnsEMBL::Transcript argument not [' . ref($transcript) .'].');
   }

   my $features = $self->fetch_all_by_ensembl_object_type('Transcript', $transcript->dbID());

   return $features;

}

=head2 fetch_all_by_Translation

  Arg [1]    : Bio::EnsEMBL::Translation
  Example    : @features = @{$regulatory_feature_adaptor->
                      fetch_all_by_Translation($Translation)};
  Description: Gets all the regulatory features associated with a
               particular Translation. Each feature only appears once.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature
  Exceptions : If arg is not of correct type.
  Caller     : ?

=cut

sub fetch_all_by_translation {

   my ($self, $translation) = @_;

   if(!ref($translation) || !$translation->isa('Bio::EnsEMBL::Translation')) {
      throw('Expected Bio::EnsEMBL::translation argument not [' . ref($translation) .'].');
    }

   my $features = $self->fetch_all_by_ensembl_object_type('Translation', $translation->dbID());

   return $features;

}

=head2 fetch_all_by_gene

  Arg [1]    : Bio::EnsEMBL::Gene
  Arg [2]    : if set, return regulatory features associated with the 
               transcripts of the gene as well.
  Example    : @features = @{$regulatory_feature_adaptor->
                      fetch_all_by_gene($gene, 1)};
  Description: Gets all the regulatory features associated with a
               particular gene, and (optionally) its transcripts.
               Each feature only appears once.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature
  Exceptions : If arg is not of correct type.
  Caller     : ?

=cut

sub fetch_all_by_gene {

   my ($self, $gene, $recursive) = @_;

   if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
      throw('Expected Bio::EnsEMBL::Gene argument not [' . ref($gene) .'].');
    }

   my @features = @{$self->fetch_all_by_ensembl_object_type('Gene', $gene->dbID())};

   # optionally add transcripts' features as well
   if ($recursive) {

     foreach my $transcript (@{$gene->get_all_Transcripts()}) {
       push @features, @{$self->fetch_all_by_transcript($transcript)};
     }

   }


   return \@features;

}

1;



#
# Ensembl module for Bio::EnsEMBL::DBSQL::AttributeAdaptor
#
# Copyright (c) 2003 EnsEMBL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AttributeAdaptor - Provides database interaction for
Bio::EnsEMBL::Attribute objects.


=head1 SYNOPSIS

  #$db is a Bio::EnsEMBL::DBSQL::DBAdaptor object:
  $attribute_adaptor = $db->get_AttributeAdaptor();

  $attributes = $attribute_adaptor->fetch_all_by_MiscFeature( $feature );

  $attributes = $attribute_adaptor->fetch_all_by_Slice( $slice );

  $attribute_adaptor->store_on_Slice( $slice, \@attributes );

  $attribute_adaptor->store_on_MiscFeature( $misc_feature, \@attributes )


=head1 DESCRIPTION


=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::AttributeAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Attribute;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Arg [...]  : Superclass args.  See Bio::EnsEMBL::DBSQL::BaseAdaptor
  Description: Instantiates a Bio::EnsEMBL::DBSQL::AttributeAdaptor
  Returntype : Bio::EnsEMBL::AttributeAdaptor
  Exceptions : none
  Caller     : DBAdaptor

=cut


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);


  # cache creation could go here
  return $self;
}



=head2 fetch_all_by_MiscFeature

  Arg [1]    : Bio::EnsEMBL::MiscFeature $mf
  Example    : @attributes = @{$attrib_adaptor->fetch_all_by_MiscFeature($mf)};
  Description: Fetches all attributes for a given MiscFeature
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : throw if incorrect arguments
               throw if provided MiscFeature does not have a dbID
  Caller     : MiscFeature

=cut

sub fetch_all_by_MiscFeature {
  my $self = shift;
  my $mf   = shift;

  if(!ref($mf) || !$mf->isa('Bio::EnsEMBL::MiscFeature')) {
    throw('MiscFeature argument is required.');
  }

  my $mfid = $mf->dbID();

  if(!defined($mfid)) {
    throw("MiscFeature must have dbID.");
  }

  my $sth = $self->prepare("SELECT at.code, at.name, at.description, " .
                           "       ma.value " .
                           "FROM misc_attrib ma, attrib_type at " .
                           "WHERE ma.misc_feature_id = ? " .
                           "AND   at.attrib_type_id = ma.attrib_type_id");

  $sth->execute($mfid);

  my $results = $self->_obj_from_sth($sth);

  $sth->finish();

  return $results;
}


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : @attributes = @{$attrib_adaptor->fetch_all_by_Slice($slice)};
  Description: Fetches all attributes for a given sequence region (which the
               passed in slice is on)
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : throw if incorrect arguments
               throw if cannot get seq_region_id from provided Slice
  Caller     : Slice

=cut

sub fetch_all_by_Slice {
  my $self = shift;
  my $slice = shift;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Slice argument is required.');
  }

  my $seq_region_id = $slice->get_seq_region_id();


  if(!defined($seq_region_id)) {
    throw("Could not get seq_region_id for provided slice: ".$slice->name());
  }

  my $sth = $self->prepare("SELECT at.code, at.name, at.description, " .
                           "       sra.value " .
                           "FROM seq_region_attrib sra, attrib_type at " .
                           "WHERE sra.seq_region_id = ? " .
                           "AND   at.attrib_type_id = sra.attrib_type_id");

  $sth->execute($seq_region_id);

  my $results = $self->_obj_from_sth($sth);

  $sth->finish();

  return $results;
}



=head2 store_on_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : listref of Bio::EnsEMBL::Attribute objects $attribs
  Example    : $attribute_adaptor->store_on_Slice($slice, \@attribs);
  Description: Stores a set of attributes on a sequence region given a
               Slice object which is on the seq_region for which attributes are
               being stored.
  Returntype : none
  Exceptions : throw if $slice argument not provided
  Caller     : general

=cut

sub store_on_Slice {
  my $self     = shift;
  my $slice    = shift;
  my $attribs = shift;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw("Slice argument expected.");
  }

  if(ref($attribs) ne 'ARRAY') {
    throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
          "argument expected.");
  }

  my $seq_region_id = $slice->get_seq_region_id();

  if(!$seq_region_id) {
    throw("Could not get seq_region_id for provided slice: ".$slice->name());
  }

  my $sth = $self->prepare( "INSERT into seq_region_attrib ".
                            "SET seq_region_id = ?, attrib_type_id = ?, ".
                            "value = ? " );

  foreach my $at ( @$attribs ) {
    if(!ref($at) && $at->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }
    my $atid = $self->_store_type( $at );
    $sth->execute( $seq_region_id, $atid, $at->value() );
  }

  return;
}




=head2 store_on_MiscFeature

  Arg [1]    : Bio::EnsEMBL::MiscFeature $feature
  Example    : $attribute_adaptor->store_on_MiscFeature($my_misc_feature,
                                                        $attributes)
  Description: Stores all the attributes on the misc feature. 
               Will duplicate things if called twice.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if provided feature is not stored in this database
  Caller     : general, MiscFeatureAdaptor

=cut

sub store_on_MiscFeature {
  my $self       = shift;
  my $feature    = shift;
  my $attributes = shift;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::MiscFeature')) {
    throw("MiscFeature argument expected");
  }

  if(ref($attributes) ne 'ARRAY') {
    throw("Reference to list of Bio::EnsEMBL::Attribute objects argument " .
          "expected");
  }

  my $db = $self->db();
  if(!$feature->is_stored($db)) {
    throw("MiscFeature is not stored in this DB - cannot store attributes.");
  }

  my $feature_id = $feature->dbID();

  my $sth = $self->prepare( "INSERT into misc_attrib ".
			    "SET misc_feature_id = ?, attrib_type_id = ?, ".
			    "value = ? " );

  for my $attrib ( @$attributes ) {
    if(!ref($attrib) && $attrib->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }
    my $atid = $self->_store_type( $attrib );
    $sth->execute( $feature_id, $atid, $attrib->value() );
  }

  return;
}




=head2 remove_from_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The Slice to remove attributes from
  Example    : $attribute_adaptor->remove_from_Slice($slice);
  Description: Removes all attributes which are stored in the database on
               the provided Slice.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if cannot obtain seq_region_id from provided Slice
  Caller     : general

=cut

sub remove_from_Slice {
  my $self = shift;
  my $slice = shift;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw("Bio::EnsEMBL::Slice argument expected.");
  }

  my $srid = $slice->get_seq_region_id();

  if(!$srid) {
    throw("Could not get seq_region_id for provided slice: ".$slice->name());
  }

  my $sth = $self->prepare("DELETE FROM seq_region_attrib " .
                           "WHERE seq_region_id = ?");

  $sth->execute($srid);

  $sth->finish();

  return;
}


=head2 remove_from_MiscFeature

  Arg [1]    : Bio::EnsEMBL::MiscFeature $mf
               The MiscFeature to remove attributes from
  Example    : $attribute_adaptor->remove_from_MiscFeature($mf);
  Description: Removes all attributes which are stored in the database on
               the provided MiscFeature.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if MiscFeature is not stored in this database
  Caller     : general

=cut

sub remove_from_MiscFeature {
  my $self = shift;
  my $mf   = shift;

  if(!ref($mf) || !$mf->isa('Bio::EnsEMBL::MiscFeature')) {
    throw("Bio::EnsEMBL::MiscFeature argument is required.");
  }

  my $db = $self->db();

  if(!$mf->is_stored($db)) {
    throw("MiscFeature is not stored in this database.");
  }

  my $sth = $db->prepare("DELETE FROM misc_attrib " .
                         "WHERE misc_feature_id = ?");

  $sth->execute($mf->dbID());

  $sth->finish();

  return;
}



# _store_type

sub _store_type {
  my $self = shift;
  my $attrib = shift;

  my $sth1 = $self->prepare
    ("INSERT IGNORE INTO attrib_type set code = ?, name = ?, ".
     "description = ?" );

  my $rows_inserted =
    $sth1->execute($attrib->code(), $attrib->name(), $attrib->desc() );

  my $atid = $sth1->{'mysql_insertid'};

  if($rows_inserted == 0) {
    # the insert failed because the code is already stored
    my $sth2 = $self->prepare
      ("SELECT attrib_type_id FROM attrib_type " .
       "WHERE code = ?");
    $sth2->execute($attrib->code());
    ($atid) = $sth2->fetchrow_array();

    $sth2->finish();

    if(!$atid) {
      throw("Could not store or fetch attrib_type code [".$attrib->code."]\n" .
	    "Wrong database user/permissions?");
    }
  }

  $sth1->finish();


  return $atid;
}




sub _obj_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($code, $name, $desc, $value);
  $sth->bind_columns(\$code, \$name, \$desc, \$value);

  my @results;
  while($sth->fetch()) {
    push @results, Bio::EnsEMBL::Attribute->new(-CODE => $code,
                                                -NAME => $name,
                                                -DESC => $desc,
                                                -VALUE => $value);
  }

  return \@results;
}


#
# Called during db destruction to clean up internal cache structures
# that result in circular references
#
sub deleteObj {
  my $self = shift;

  #break circular db <-> adaptor references
  $self->SUPER::deleteObj();

  #break circular object <-> adaptor references

}


1;

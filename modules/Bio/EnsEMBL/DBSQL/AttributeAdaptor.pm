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

  $attribute_types = $attribute_adaptor->fetch_all_types();

  $attribute_adaptor->store_on_Slice( $slice, $attributes );

  $attribute_adaptor->store_on_MiscFeature( $misc_feature, $attributes )


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

  Arg [1]    : none, string, int, Bio::EnsEMBL::Example $formal_parameter_name
    Additional description lines
    list, listref, hashref
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut



=head2 fetch_all_by_Slice

  Arg [1]    : none, string, int, Bio::EnsEMBL::Example $formal_parameter_name
    Additional description lines
    list, listref, hashref
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut



=head2 fetch_all_types

  Arg [1]    : none, string, int, Bio::EnsEMBL::Example $formal_parameter_name
    Additional description lines
    list, listref, hashref
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut

sub fetch_all_types {
  my $self = shift;

  return $all_types;
}






=head2 store_on_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut


sub store_on_Slice {
  my $self = shift;
  my $slice = shift;

  # need to access directly to prevent lazy loading
  my $attributes = $slice->{'???????'};

  if( ! defined $attributes ) {
    return;
  }

  my $seq_region_id = $slice->get_seq_region_id();

  if(  ! defined $seq_region_id ) {
    throw( );
  }

  my $sth = $self->prepare( "INSERT into seq_region_attrib ".
			    "SET seq_region_id = ?, attrib_type_id = ?, ".
			    "value = ? " );

  for my $attrib ( @$attributes ) {
    my $attrib_id = $self->_store_type( $attrib );
    $sth->excute( $seq_region_id, $attrib_id, $attrib->value() );
  }

  return;

}




=head2 store_on_MiscFeature

  Arg [1]    : Bio::EnsEMBL::MiscFeature $feature
  Example    : $attribute_adaptor->store_on_MiscFeature( $my_misc_feature )
  Description: Stores all the attributes on the misc feature. 
               Will duplicate things if called twice.
  Returntype : none
  Exceptions : on database problems.
  Caller     : general, MiscFeatureAdaptor

=cut

sub store_on_MiscFeature {
  my $self = shift;
  my $feature = shift;

  # need to access directly to prevent lazy loading
  my $attributes = $feature->{'attributes'};

  if( ! defined $attributes ) {
    return;
  }

  my $feature_id = $feature->dbID();

  if(  ! defined $feature_id ) {
    throw("dbID in feature is not set.".
	  "Need to store feature before storeing attributes.\n" );
  }

  my $sth = $self->prepare( "INSERT into misc_attrib ".
			    "SET misc_feature_id = ?, attrib_type_id = ?, ".
			    "value = ? " );

  for my $attrib ( @$attributes ) {
    my $attrib_id = $self->_store_type( $attrib );
    $sth->excute( $feature_id, $attrib_id, $attrib->value() );
  }

  return;
}



# _store_type 

sub _store_type {
  my $self = shift;
  my $attrib = shift;

  
  my $sth = $self->prepare
    ("INSERT IGNORE INTO attrib_type set code = ?, name = ?, ".
     "description = ?" );

  $sth->execute($attrib->code(), $attrib->name(), $attrib->description() );

  my $atid = $sth->{'mysql_insertid'};

  $sth->finish();

  if( $self->db->db_handle->{'mysql_info'} == 0 ){
    # the insert failed because the code is already stored
    $sth = $self->prepare
      ("SELECT attrib_type_id FROM attrib_type " .
       "WHERE code = ?");
    $sth->execute($code);
    ($atid) = $sth->fetchrow_array();

    $sth->finish();

    if(!$atid) {
      throw("Could not store or fetch attrib_type code [$code]\n" .
	    "Wrong permissions?");
    }
  }

  return $atid;
}


sub _obj_from_sth {
  my $self = shift;
  my $sth = shift;

  return $results;
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

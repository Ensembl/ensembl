#
# Ensembl module for Bio::EnsEMBL::MiscFeature
#
# Copyright (c) 2003 EnsEMBL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::MiscFeature - A miscelaneous feature with arbitrary features and
associations.

=head1 SYNOPSIS

  use Bio::EnsEMBL::MiscFeature;
  use Bio::EnsEMBL::MiscSet;
  use Bio::EnsEMBL::Attribute;

  my $mfeat = Bio::EnsEMBL::MiscFeature->new(-START  => 1200,
                                             -END    => 100_000,
                                             -STRAND => 1,
                                             -SLICE  => $slice);

  #  Can add attributes to the misc feature and associate with various sets
  my $clone_set = Bio::EnsEMBL::MiscSet->new(-CODE => 'clone',
                                         -NAME => '1MB clone set',
                                         -DESCRIPTION => '1MB CloneSet');


  my $tiling_path_set = Bio::EnsEMBL::MiscSet->new(-CODE => 'tilingpath',
                                                   -NAME => 'tiling path set');


  my $attrib1 = Bio::EnsEMBL::Attribute->new(-VALUE => 'RLX12451',
                                            -CODE => 'name',
                                            -NAME => 'name');

  my $attrib2 = Bio::EnsEMBL::Attribute->new(-VALUE => '4',
                                            -CODE  => 'version',
                                           -NAME  => 'version');

  my $attrib3 = Bio::EnsEMBL::Attribute->new(-VALUE => 'AL42131.4',
                                            -CODE  => 'synonym',
                                            -NAME  => 'synonym');

  # can associate a misc feature with any number of sets

  $mfeat->add_MiscSet($clone_set);
  $mfeat->add_MiscSet($tiling_path_set);

  # can add arbitrary attributes to a misc feature

  $mfeat->add_Attribute($attrib1);
  $mfeat->add_Attribute($attrib2);
  $mfeat->add_Attribute($attrib3);

  my ($name_attrib)   = @{$mfeat->get_all_Attributes('name')};
  my @all_attribs     = @{$mfeat->get_all_Attributes()};

  my @all_sets = @{$mfeat->get_all_MiscSets()};
  my ($clone_set) = @{$mfeat->get_all_CloneSets('clone')};


  # Can do normal feature operations as well
  $mfeat = $mfeat->transform('supercontig');
  print $mfeat->slice->seq_region_name, ' ', $mfeat->start, '-', $mfeat->end;


=head1 DESCRIPTION

MiscFeatures are extremely general  features with a location, and an
arbitrary group of attributes.  They are grouped with other features of the
same 'type' through the use of MiscSets (see Bio::EnsEMBL::MiscSet).
Attributes are attached in the fom of Bio::EnsEMBL::Attribute objects.
See Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor for ways to fetch or store
MiscFeatures.

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=cut


use strict;
use warnings;

package Bio::EnsEMBL::MiscFeature;

use Bio::EnsEMBL::Feature;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Feature);



# new is inherited from superclass

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless($hashref, $class);
}



=head2 add_Attribute

  Arg [1]    : Bio::EnsEMBL::Attribute $attribute
  Example    : $misc_feature->add_attribute($attribute);
  Description: Adds an attribute to this misc. feature
  Returntype : none
  Exceptions : throw on wrong argument type
  Caller     : general

=cut

sub add_Attribute {
  my ($self, $attrib) = @_;

  if( ! defined $attrib || ! $attrib->isa( "Bio::EnsEMBL::Attribute" )) {
    throw( "You have to provide a Bio::EnsEMBL::Attribute, not a [$attrib]" );
  }

  $self->{'attributes'} ||= [];
  push @{$self->{'attributes'}}, $attrib
}



=head2 add_MiscSet

  Arg [1]    : Bio::EnsEMBL::MiscSet $set
               The set to add
  Example    : $misc_feature->add_MiscSet(Bio::EnsEMBL::MiscSet->new(...));
  Description: Associates this MiscFeature with a given Set.
  Returntype : none
  Exceptions : throw if the set arg is not provided,
               throw if the set to be added does not have a code
  Caller     : general

=cut


sub add_MiscSet {
  my $self = shift;
  my $miscSet = shift;

  if(!$miscSet || !ref($miscSet) || !$miscSet->isa('Bio::EnsEMBL::MiscSet')) {
    throw('Set argument must be a Bio::EnsEMBL::MiscSet');
  }

  $self->{'miscSets'} ||= [];

  push( @{$self->{'miscSets'}}, $miscSet );
}



=head2 get_all_MiscSets

  Arg [1]    : optional string $code
               The code of the set to retrieve
  Example    : $set = $misc_feature->get_all_MiscSets($code);
  Description: Retrieves a set that this feature is associated with via its
               code. Can return empty lists. Usually returns about one elements lists.
  Returntype : listref of Bio::EnsEMBL::MiscSet
  Exceptions : throw if the code arg is not provided
  Caller     : general

=cut


sub get_all_MiscSets {
  my $self = shift;
  my $code = shift;

  $self->{'miscSets'} ||= [];
  if( defined $code ) {
    my @results = grep { uc($_->code())eq uc( $code ) } @{$self->{'miscSets'}};
    return \@results;
  } else {
    return $self->{'miscSets'};
  }
}


=head2 get_all_Attributes

  Arg [1]    : optional string $code
               The code of the Attribute objects to retrieve
  Example    : @attributes = $misc_feature->get_all_Attributes('name');
  Description: Retrieves a list of Attribute objects for given code or all
               of the associated Attributes.
  Returntype : listref of Bio::EnsEMBL::Attribute
  Exceptions : 
  Caller     : general

=cut

sub get_all_Attributes {
  my $self = shift;
  my $code = shift;

  my @results;
  my $result;

  if( defined $code ) {
    @results = grep { uc( $_->code() ) eq uc( $code )} @{$self->{'attributes'}};
    return \@results;
  } else {
    return $self->{'attributes'};
  }
}

sub get_all_attribute_values {
  my $self = shift;
  my $code = shift;
  my @results = map { uc( $_->code() ) eq uc( $code ) ? $_->value : () } @{$self->{'attributes'}};
  return \@results;
}

sub get_scalar_attribute {
  my $self = shift;
  my $code = shift;
  my @results;
  if( $code eq 'name' ) {
    @results = grep { uc( $_->code() ) eq 'NAME' || uc( $_->code() ) eq 'NON_REF' } @{$self->{'attributes'}};
  } else {
    @results = grep { uc( $_->code() ) eq uc( $code )} @{$self->{'attributes'}};
  }
  return @results ? $results[0]->value() : '';
}

=head2 display_id

  Arg [1]    : none
  Example    : print $kb->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For misc_features this is the first
               name or synonym attribute or '' if neither are defined.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code

=cut

sub display_id {
  my $self = shift;
  my ($attrib) = @{$self->get_all_Attributes('name')};
  ($attrib) =  @{$self->get_all_Attributes('non_ref')} if(!$attrib);
  ($attrib) =  @{$self->get_all_Attributes('synonym')} if(!$attrib);
  if( defined $attrib ) {
    return $attrib->value();
  } else {
    return '';
  }
}




1;

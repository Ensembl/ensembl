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

  my $mfeat = Bio::EnsEMBL::MiscFeature->new(-START  => 1200,
                                             -END    => 100_000,
                                             -STRAND => 1,
                                             -SLICE  => $slice);

  #  Can add attributes to the misc feature and associate with various sets
  $clone_set->code('clone');
  $tiling_path_set->core('tilingpath');

  $mfeat->add_attribute('name', 'RLX12451');
  $mfeat->add_attribute('version', '4');
  $mfeat->add_attribute('synonym', 'AL42131.4');
  $mfeat->add_attrubute('synonym', 'Z199311');
  $mfeat->add_set($clone_set);
  $mfeat->add_set($tiling_path_set);

  my ($name)   = $mfeat->get_attribute('name');
  my @synonyms = $mfeat->get_attribute('synonym');

  my @attrib_types = $mfeat->get_attribute_types();
  foreach my $attrib_type (@attrib_types) {
    print join(',', $mfeat->get_attribute($attrib_type));
  }

  $clone_set = $mfeat->get_set('clone');
  $tiling_path_set = $mfeat->get_set('tiling_path');

  my @set_codes = $mfeat->get_set_codes();
  foreach my $set_code (@set_codes) {
    print $mfeat->get_set->name(), "\n";
  }


  # Can do normal feature operations as well
  $mfeat = $mfeat->transform('supercontig');
  print $mfeat->slice->seq_region_name, ' ', $mfeat->start, '-', $mfeat->end;


=head1 DESCRIPTION

MiscFeatures are extremely general  features with a location, and an
arbitrary group of attributes.  They are grouped with other features of the
same 'type' through the use of MiscSets (see Bio::EnsEMBL::MiscSet).

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=cut


use strict;
use warnings;

use Bio::EnsEMBL::Feature;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use vars qw(@ISA);


# new is inherited from superclass

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless($hashref), $class;
}




=head2 add_attribute

  Arg [1]    : string $type
               The type of attribute to add
  Arg [2]    : string $value
               The value of the attribute to add
  Example    : $misc_feature->add_attribute('synonym', 'AL124141.1');
  Description: Adds an attribute of a given type to this misc. feature
  Returntype : none
  Exceptions : throw if type arg is not provided
  Caller     : general

=cut

sub add_attribute {
  my ($self, $type, $value) = @_;

  throw('Type argument is required') if(!$type);

  $self->{'attributes'}->{$type} ||= [];
  push @{$self->{'attributes'}->{$type}}, $value
}


=head2 add_set

  Arg [1]    : Bio::EnsEMBL::MiscSet $set
               The set to add
  Example    : $misc_feature->add_set(Bio::EnsEMBL::MiscSet->new(...));
  Description: Associates this MiscFeature with a given Set.
  Returntype : none
  Exceptions : throw if the set arg is not provided,
               throw if the set to be added does not have a code
  Caller     : general

=cut

sub add_set {
  my ($self, $set) = @_;

  if(!$set || !ref($set) || !$set->isa('Bio::EnsEMBL::MiscSet')) {
    throw('Set argument must be a Bio::EnsEMBL::MiscSet');
  }

  my $code = $set->code();

  throw('Set must be associated with a code to be added') if(!$code);

  $self->{'sets'}->{$code} = $set;
}



=head2 get_attribute_types

  Arg [1]    : none
  Example    : @attrib_types = $misc_feature->get_attribute_types();
  Description: returns a list of attribute types that this feature has
  Returntype : list of strings
  Exceptions : none
  Caller     : general

=cut

sub get_attribute_types {
  my $self = shift;
  $self->{'attributes'} ||= {};

  return keys %{$self->{'attributes'}};
}

=head2 get_set_codes

  Arg [1]    : none
  Example    : @set_codes = $misc_feature->get_set_codes();
  Description: returns a list of codes for the sets this feature is associated
               with
  Returntype : list of strings
  Exceptions : none
  Caller     : general

=cut

sub get_set_codes {
  my $self = shift;
  $self->{'sets'} ||= {};
  return keys %{$self->{'sets'}};
}


=head2 get_set

  Arg [1]    : string $code
               The code of the set to retrieve
  Example    : $set = $misc_feature->get_set($code);
  Description: Retrieves a set that this feature is associated with via its
               code.  Undef is returned if the feature is not associated with
               the requested set.
  Returntype : Bio::EnsEMBL::MiscSet or undef
  Exceptions : throw if the code arg is not provided
  Caller     : general

=cut

sub get_set {
  my $self = shift;
  my $code = shift;

  throw('Code arg is required.') if (!$code);

  return $self->{'sets'}->{$code};
}



=head2 get_attribute

  Arg [1]    : string $type
               The type of the attribute values to retrieve
  Example    : @values = $misc_feature->get_attribute('name');
  Description: Retrieves a list of values for a particular attribute type
  Returntype : list of strings
  Exceptions : thrown if a type arg is not provided
  Caller     : general

=cut

sub get_attribute {
  my $self = shift;
  my $type = shift;

  throw('Type arg is required.') if(!$type);

  return $self->{'attributes'}->{$type} || ();
}

1;

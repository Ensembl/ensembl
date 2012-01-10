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

Bio::EnsEMBL::MiscSet - This is a set representing a classification of
a group of miscellaneuos features.

=head1 SYNOPSIS

  use Bio::EnsEMBL::MiscSet;

  my $misc_set = Bio::EnsEMBL::MiscSet->new(
    1234, $adaptor, 'tilepath',
    'Assembly Tiling Path',
    'The tiling path of clones', 1e6
  );

  my $misc_feature->add_set($misc_set);

=head1 DESCRIPTION

MiscSets represent classsifications or groupings of MiscFeatures.
Features are classified into sets essentially to define what they are
and how they may be used.  Generally MiscFeatures are retrieved on
the basis of their associated sets. See Bio::EnsEMBL::MiscFeature,
Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor.

Note that MiscSets and MiscFeatures were formerly known as MapSets and
MapFrags

=head1 METHODS

=cut

package Bio::EnsEMBL::MiscSet;

use strict;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [1]    : int $misc_set_id
               The internal identifier for this misc set
  Arg [2]    : string $code
               The unique code which identifies this set type
  Arg [3]    : string $name
               The human readable name of this set
  Arg [4]    : string $desc
               The description of this set
  Arg [5]    : int $max_len
               The maximum length of features of this mapset
  Example    : $set = new Bio::EnsEMBL::MiscSet(1234, 'tilepath',
                                                'Assembly Tiling Path',
                                                'The tiling path of clones',
                                                1e6);
  Description: Instantiates a Bio::EnsEMBL::MiscSet
  Returntype : Bio::EnsEMBL::MiscSet
  Exceptions : none
  Caller     : MiscFeatureAdaptor
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my($code, $name, $desc, $max_len) =
    rearrange([qw(CODE NAME DESCRIPTION LONGEST_FEATURE)], @_);

  $self->{'code'} = $code;
  $self->{'name'} = $name;
  $self->{'description'} = $desc;
  $self->{'longest_feature'} = $max_len;

  return $self;
}

=head2 code

  Arg [1]    : string $newval (optional) 
               The new value to set the code attribute to
  Example    : $code = $obj->code()
  Description: Getter/Setter for the code attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub code{
  my $self = shift;
  $self->{'code'} = shift if(@_);
  return $self->{'code'};
}


=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the code attribute to
  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}


=head2 description

  Arg [1]    : string $newval (optional)
               The new value to set the description attribute to
  Example    : $description = $obj->description()
  Description: Getter/Setter for the description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description{
  my $self = shift;
  $self->{'description'} = shift if(@_);
  return $self->{'description'};
}


=head2 longest_feature

  Arg [1]    : int $newval (optional) 
               The new value to set the longest_feature attribute to
  Example    : $longest_feature = $obj->longest_feature()
  Description: Getter/Setter for the longest_feature attribute
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub longest_feature{
  my $self = shift;
  $self->{'longest_feature'} = shift if(@_);
  return $self->{'longest_feature'};
}


1;

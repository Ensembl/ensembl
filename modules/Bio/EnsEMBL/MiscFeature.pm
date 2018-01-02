=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::MiscFeature - A miscelaneous feature with arbitrary features and
associations.

=head1 SYNOPSIS

  use Bio::EnsEMBL::MiscFeature;
  use Bio::EnsEMBL::MiscSet;
  use Bio::EnsEMBL::Attribute;

  my $mfeat = Bio::EnsEMBL::MiscFeature->new(
    -START  => 1200,
    -END    => 100_000,
    -STRAND => 1,
    -SLICE  => $slice
  );

  # Can add attributes to the misc feature and associate with various
  # sets
  my $clone_set = Bio::EnsEMBL::MiscSet->new(
    -CODE        => 'clone',
    -NAME        => '1MB clone set',
    -DESCRIPTION => '1MB CloneSet'
  );

  my $tiling_path_set = Bio::EnsEMBL::MiscSet->new(
    -CODE => 'tilingpath',
    -NAME => 'tiling path set'
  );

  my $attrib1 = Bio::EnsEMBL::Attribute->new(
    -VALUE => 'RLX12451',
    -CODE  => 'name',
    -NAME  => 'name'
  );

  my $attrib2 = Bio::EnsEMBL::Attribute->new(
    -VALUE => '4',
    -CODE  => 'version',
    -NAME  => 'version'
  );

  my $attrib3 = Bio::EnsEMBL::Attribute->new(
    -VALUE => 'AL42131.4',
    -CODE  => 'synonym',
    -NAME  => 'synonym'
  );

  # can associate a misc feature with any number of sets

  $mfeat->add_MiscSet($clone_set);
  $mfeat->add_MiscSet($tiling_path_set);

  # can add arbitrary attributes to a misc feature

  $mfeat->add_Attribute($attrib1);
  $mfeat->add_Attribute($attrib2);
  $mfeat->add_Attribute($attrib3);

  my ($name_attrib) = @{ $mfeat->get_all_Attributes('name') };
  my @all_attribs = @{ $mfeat->get_all_Attributes() };

  my @all_sets = @{ $mfeat->get_all_MiscSets() };
  my ($clone_set) = @{ $mfeat->get_all_CloneSets('clone') };


  # Can do normal feature operations as well
  $mfeat = $mfeat->transform('supercontig');
  print $mfeat->slice->seq_region_name, ' ', $mfeat->start, '-',
    $mfeat->end;


=head1 DESCRIPTION

MiscFeatures are extremely general features with a location, and an
arbitrary group of attributes.  They are grouped with other features of
the same 'type' through the use of MiscSets (see Bio::EnsEMBL::MiscSet).
Attributes are attached in the fom of Bio::EnsEMBL::Attribute objects.
See Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor for ways to fetch or store
MiscFeatures.

=cut


package Bio::EnsEMBL::MiscFeature;

use strict;
use warnings;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Feature);

=head2 new

  Arg [-SLICE]: Bio::EnsEMBL::SLice - Represents the sequence that this
                feature is on. The coordinates of the created feature are
                relative to the start of the slice.
  Arg [-START]: The start coordinate of this feature relative to the start
                of the slice it is sitting on.  Coordinates start at 1 and
                are inclusive.
  Arg [-END]  : The end coordinate of this feature relative to the start of
                the slice it is sitting on.  Coordinates start at 1 and are
                inclusive.
  Arg [-STRAND]: The orientation of this feature.  Valid values are 1,-1,0.
  Arg [-SEQNAME] : A seqname to be used instead of the default name of the
                of the slice.  Useful for features that do not have an
                attached slice such as protein features.
  Arg [-dbID]   : (optional) internal database id
  Arg [-ADAPTOR]: (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor
  Example    : $feature = Bio::EnsEMBL::MiscFeature->new(-start    => 1,
                                                     -end      => 100,
                                                     -strand   => 1,
                                                     -slice    => $slice,
                                                     -analysis => $analysis);
  Description: Constructs a new Bio::EnsEMBL::Feature.  Generally subclasses
               of this method are instantiated, rather than this class itself.
  Returntype : Bio::EnsEMBL::MiscFeature
  Exceptions : Thrown on invalid -SLICE, -ANALYSIS, -STRAND ,-ADAPTOR arguments
  Caller     : general, subclass constructors
  Status     : Stable

=cut


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->{'attributes'} = [];

  return $self;
}



=head2 add_Attribute

  Arg [1]    : Bio::EnsEMBL::Attribute $attribute
  Example    : $misc_feature->add_attribute($attribute);
  Description: Adds an attribute to this misc. feature
  Returntype : none
  Exceptions : throw on wrong argument type
  Caller     : general
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Example    : @attributes = @{ $misc_feature->get_all_Attributes('name') };
  Description: Retrieves a list of Attribute objects for given code or all
               of the associated Attributes.
  Returntype : listref of Bio::EnsEMBL::Attribute
  Exceptions : 
  Caller     : general
  Status     : Stable

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

=head2 get_all_attribute_values

  Arg [1]    : string $code
               The code of the Attribute object values to retrieve
  Example    : @attributes_vals = @{$misc_feature->get_all_attribute_values('name')};
  Description: Retrieves a list of Attribute object values for given code or all
               of the associated Attributes.
  Returntype : listref of values
  Exceptions : 
  Caller     : general
  Status     : Stable

=cut

sub get_all_attribute_values {
  my $self = shift;
  my $code = shift;
  my @results = map { uc( $_->code() ) eq uc( $code ) ? $_->value : () } 
                @{$self->{'attributes'}};
  return \@results;
}

=head2 get_scalar_attribute

  Arg [1]    : string $code
               The code of the Attribute object values to retrieve
  Example    : $vals = $misc_feature->get_scalar_attribute('name');
  Description: Retrieves a value for given code or all
               of the associated Attributes.
  Returntype : scalar value
  Exceptions : 
  Caller     : general
  Status     : Stable


=cut
  

sub get_scalar_attribute {
  my $self = shift;
  my $code = shift;
  my @results = grep { uc( $_->code() ) eq uc( $code )} @{$self->{'attributes'}};
  return @results ? $results[0]->value() : '';
}

sub get_first_scalar_attribute {
  my $self = shift;
  foreach my $code ( @_ ) {
    my @results = grep { uc( $_->code() ) eq uc( $code )} @{$self->{'attributes'}};
    return $results[0]->value() if @results;
  }
  return '';
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
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  my ($attrib) = @{$self->get_all_Attributes('name')};
  ($attrib) =  @{$self->get_all_Attributes('synonym')} if(!$attrib);
  if( defined $attrib ) {
    return $attrib->value();
  } else {
    return '';
  }
}

=head2 summary_as_hash

  Example    : my $hash = $misc_feature->summary_as_hash();
  Description: Generates a HashRef compatible with GFFSerializer. Adds
               all attribute key value pairs plus MiscSet codes and names
  Returntype : Hash
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub summary_as_hash {
  my ($self) = @_;
  my $hash = $self->SUPER::summary_as_hash();
  my $attributes = $self->get_all_Attributes();
  foreach my $attr (@{$attributes}) {
    $hash->{$attr->code()} = $attr->value();
  }
  my $misc_sets = $self->get_all_MiscSets();
  foreach my $set (@{$misc_sets}) {
    push(@{$hash->{misc_set_code}},$set->code());
    push(@{$hash->{misc_set_name}},$set->name());
  }
  return $hash;
}


1;

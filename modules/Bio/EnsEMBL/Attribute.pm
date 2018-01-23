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

Bio::EnsEMBL::Attribute - A generic Attribute class.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Attribute;

  my $attribute = Bio::EnsEMBL::Attribute->new
       (-CODE => 'myCode',
        -NAME => 'My Attribute',
        -DESCRIPTION => 'This is my attribute description.',
        -VALUE => '10023');

  print $attrib->name(), "\n";
  print $attrib->code(), "\n";
  print $attrib->description(), "\n";
  print $attrib->value(), "\n";

=head1 DESCRIPTION

This is a generic attribute class used to represent attributes
associated with seq_regions (and their Slices) and MiscFeatures.

=head1 SEE ALSO

Bio::EnsEMBL::Slice
Bio::EnsEMBL::MiscFeature
Bio::EnsEMBL::DBSQL::AttributeAdaptor

=cut

package Bio::EnsEMBL::Attribute;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

=head2 new

  Arg [-CODE]        : string - the code for this attribute
  Arg [-NAME]        : string - a human readable name for this attribute
  Arg [-DESCRIPTION] : string - a description for this attribute
  Arg [-VALUE]       : value  - the value of this attribute
  Example            :   my $attribute = Bio::EnsEMBL::Attribute->new
                         (-CODE => 'myCode',
                          -NAME => 'My Attribute',
                          -DESCRIPTION => 'This is my attribute description.',
                          -VALUE => '10023');
  Description        : Constructor.  Instantiates a Bio::EnsEMBL::Attribute object.
  Returntype         : Bio::EnsEMBL::Attribute
  Exceptions         : none
  Caller             : general
  Status             : Stable

=cut


sub new {
  my $caller = shift;

  # allow to be called as class or object method
  my $class = ref($caller) || $caller;

  my ($code, $name, $desc, $value) =
    rearrange([qw(CODE NAME DESCRIPTION VALUE)], @_);

  return bless {'code'    => $code,
                'name'    => $name,
                'description' => $desc,
                'value'   => $value}, $class;
}

=head2 new_fast

  Arg [1]    : hashref to be blessed
  Description: Construct a new Bio::EnsEMBL::Attribute using the hashref.
  Exceptions : none
  Returntype : Bio::EnsEMBL::Attribute
  Caller     : general, subclass constructors
  Status     : Stable

=cut


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  my $self = bless $hashref, $class;
  weaken($self->{adaptor})  if ( ! isweak($self->{adaptor}) );
  return $self;
}


=head2 code

  Arg [1]    : string $code (optional)
  Example    : $code = $attribute->code();
  Description: Getter/Setter for code attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub code {
  my $self = shift;
  $self->{'code'} = shift if(@_);
  return $self->{'code'};
}


=head2 name

  Arg [1]    : string $name (optional)
  Example    : $name = $attribute->name();
  Description: Getter/Setter for name attribute
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

  Arg [1]    : string $description (optional)
  Example    : $description = $attribute->description();
  Description: Getter/Setter for description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
  my $self = shift;
  $self->{'description'} = shift if(@_);
  return $self->{'description'};
}


=head2 value

  Arg [1]    : string $value (optional)
  Example    : $value = $attribute->value();
  Description: Getter/Setter for value attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub value {
  my $self = shift;
  $self->{'value'} = shift if(@_);
  return $self->{'value'};
}


1;

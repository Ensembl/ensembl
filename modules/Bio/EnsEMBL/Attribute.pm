#
# Ensembl module for Bio::EnsEMBL::Attribute
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

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

Also see B<Bio::EnsEMBL::Slice>, B<Bio::EnsEMBL::MiscFeature> and
B<Bio::EnsEMBL::DBSQL::AttributeAdaptor>.

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Attribute;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);

=head2 new

  Arg [-CODE]   : string - the code for this attribute
  Arg [-NAME]   : string - a human readable name for this attribute
  Arg [-DESCRIPTION]   : string - a description for this attribute
  Arg [-VALUE]  : value  - the value of this attribute
  Example    :   my $attribute = Bio::EnsEMBL::Attribute->new
                      (-CODE => 'myCode',
                       -NAME => 'My Attribute',
                       -DESCRIPTION => 'This is my attribute description.',
                       -VALUE => '10023');
  Description: Constructor.  Instantiates a Bio::EnsEMBL::Attribute object.
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub new {
  my $caller = shift;

  # allow to be called as class or object method
  my $class = ref($caller) || $caller;

  my ($code, $name, $desc, $value) =
    rearrange([qw(CODE NAME DESCRIPTION VALUE)], @_);

  return bless {'code'    => $code,
                'name'    => $name,
                'description'    => $desc,
                'value'   => $value}, $class;
}


=head2 code

  Arg [1]    : string $code (optional)x
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

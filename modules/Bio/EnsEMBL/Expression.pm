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

Bio::EnsEMBL::Expression - A generic Expression class.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Expression;

  my $expression = Bio::EnsEMBL::Expression->new
       (-NAME => 'My Tissue',
        -DESCRIPTION => 'This is my tissue description.',
        -VALUE => '0.8');

  print $expression->name(), "\n";
  print $expression->description(), "\n";
  print $expression->value(), "\n";

=head1 DESCRIPTION

This is a generic attribute class used to represent attributes
associated with seq_regions (and their Slices) and MiscFeatures.

=head1 SEE ALSO

Bio::EnsEMBL::DBSQL::ExpressionAdaptor

=cut

package Bio::EnsEMBL::Expression;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

=head2 new

  Arg [-NAME]        : string - the name for this tissue
  Arg [-DESCRIPTION] : string - a description for this tissue
  Arg [-VALUE]       : value  - the expression value for the tissue in a given object
  Example            :   my $expression = Bio::EnsEMBL::Expression->new
                           (-NAME => 'My Tissue',
                            -DESCRIPTION => 'This is my tissue description.',
                            -VALUE => '0.8');
  Description        : Constructor.  Instantiates a Bio::EnsEMBL::Expression object.
  Returntype         : Bio::EnsEMBL::Expression
  Exceptions         : none
  Caller             : general
  Status             : Stable

=cut


sub new {
  my $caller = shift;

  # allow to be called as class or object method
  my $class = ref($caller) || $caller;

  my ($name, $desc, $object, $value) =
    rearrange([qw(NAME DESCRIPTION OBJECT VALUE)], @_);

  return bless {'name'    => $name,
                'description' => $desc,
                'object' => $object,
                'value'   => $value}, $class;
}

=head2 new_fast

  Arg [1]    : hashref to be blessed
  Description: Construct a new Bio::EnsEMBL::Expression using the hashref.
  Exceptions : none
  Returntype : Bio::EnsEMBL::Expression
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

=head2 object

  Arg [1]    : string $object (optional)
  Example    : $object = $attribute->object();
  Description: Getter/Setter for object expression
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub object {
  my $self = shift;
  $self->{'object'} = shift if(@_);
  return $self->{'object'};
}



1;

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

Bio::EnsEMBL::IdMapping::Entry - object representing a ScoredMappingMatrix entry

=head1 SYNOPSIS

=head1 DESCRIPTION

This object represents a ScoredMappingMatrix entry. It is defined by a
pair of a source and target object's internal Id and a score for this
mapping.

=head1 METHODS

  new
  new_fast
  source
  target
  score
  to_string

=cut

package Bio::EnsEMBL::IdMapping::Entry;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 new

  Example     : my $entry = Bio::EnsEMBL::IdMapping::Entry->new();
  Description : Constructor. This is a no-argument constructor, so you need to
                populate the object manually. Rarely used since in most cases
                new_fast() is preferred.
  Return type : a Bio::EnsEMBL::IdMapping::Entry object
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = [];
  bless ($self, $class);

  return $self;
}


=head2 new_fast

  Arg[1]      : Arrayref $array_ref - the arrayref to bless into the Entry
                object 
  Example     : my $entry = Bio::EnsEMBL::IdMapping::Entry->new_fast([
                  $source_gene->id, $target_gene->id, 0.9]);
  Description : Fast constructor.
  Return type : a Bio::EnsEMBL::IdMapping::Entry object
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new_fast {
  my $class = shift;
  my $array_ref = shift;
  return bless $array_ref, $class;
}


=head2 source

  Arg[1]      : (optional) Int - source object's internal Id
  Description : Getter/setter for source object's internal Id.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub source {
  my $self = shift;
  $self->[0] = shift if (@_);
  return $self->[0];
}


=head2 target

  Arg[1]      : (optional) Int - target object's internal Id
  Description : Getter/setter for target object's internal Id.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub target {
  my $self = shift;
  $self->[1] = shift if (@_);
  return $self->[1];
}


=head2 score

  Arg[1]      : (optional) Float - a score
  Description : Getter/setter for score for the mapping between source and
                target object.
  Return type : Float
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub score {
  my $self = shift;
  $self->[2] = shift if (@_);
  return $self->[2];
}


=head2 to_string

  Example     : print LOG $entry->to_string, "\n";
  Description : Returns a string representation of the Entry object. Useful for
                debugging and logging.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub to_string {
  my $self = shift;
  return sprintf('%-10s%-10s%-5.6f', $self->source, $self->target, $self->score);
}


1;


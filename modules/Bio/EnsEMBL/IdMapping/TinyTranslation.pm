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

Bio::EnsEMBL::IdMapping::TinyTranslation - lightweight translation object

=head1 SYNOPSIS

  if ( my $tl = $tr->translation ) {
    my $lightweight_tl =
      Bio::EnsEMBL::IdMapping::TinyTranslation->new_fast( [
        $tl->dbID,          $tl->stable_id,
        $tl->version,       $tl->created_date,
        $tl->modified_date, $tr->dbID,
        $tr->translate->seq, ( $tr->is_known ? 1 : 0 ),
      ] );
  }

=head1 DESCRIPTION

This is a lightweight translation object for the stable Id mapping. See
the documentation in TinyFeature for general considerations about its
design.

=head1 METHODS

  transcript_id
  seq
  is_known

=cut

package Bio::EnsEMBL::IdMapping::TinyTranslation;

# internal data structure (array indices):
#
#  0-4 see TinyFeature
#  5  transcript_id
#  6  seq
#  7  is_known


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::TinyFeature;
our @ISA = qw(Bio::EnsEMBL::IdMapping::TinyFeature);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 transcript_id

  Arg[1]      : (optional) Int - the transcript internal Id ("dbID")
  Description : Getter/setter for the transcript internal Id this translation is
                attached to.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub transcript_id {
  my $self = shift;
  $self->[5] = shift if (@_);
  return $self->[5];
}


=head2 seq

  Arg[1]      : (optional) String - the translation's sequence
  Description : Getter/setter for the translation's sequence.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq {
  my $self = shift;
  $self->[6] = shift if (@_);
  return $self->[6];
}


=head2 is_known

  Arg[1]      : (optional) Boolean - the translation's "known" status
  Description : Getter/setter for the translation's "known" status.
  Return type : Boolean
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_known {
  my $self = shift;
  $self->[7] = shift if (@_);
  return $self->[7];
}


1;


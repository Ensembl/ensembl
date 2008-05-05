package Bio::EnsEMBL::IdMapping::TinyTranslation;

=head1 NAME

Bio::EnsEMBL::IdMapping::TinyTranslation - lightweight translation object

=head1 SYNOPSIS

if (my $tl = $tr->translation) {
  my $lightweight_tl = Bio::EnsEMBL::IdMapping::TinyTranslation->new_fast([
      $tl->dbID,
      $tl->stable_id,
      $tl->version,
      $tl->created_date,
      $tl->modified_date,
      $tr->dbID,
      $tr->translate->seq,
      ($tr->is_known ? 1 : 0),
  ]);
}

=head1 DESCRIPTION

This is a lightweight translation object for the stable Id mapping. See the
documentation in TinyFeature for general considerations about its design.

=head1 METHODS

transcript_id
seq
is_known

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


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


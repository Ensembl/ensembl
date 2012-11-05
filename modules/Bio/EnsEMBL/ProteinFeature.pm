
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

Bio::EnsEMBL::ProteinFeature

=head1 SYNOPSIS

  my $feature = Bio::EnsEMBL::ProteinFeature->new(
    -start    => $start,
    -end      => $end,
    -hstart   => $hit_start,
    -hend     => $hit_end,
    -hseqname => $hit_name
  );

=head1 DESCRIPTION

ProteinFeature objects represent domains or other features of interest
on a peptide sequence.

=head1 METHODS

=cut

package Bio::EnsEMBL::ProteinFeature;

use strict;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::FeaturePair);

=head2 new

  Arg [IDESC]           : (optional) string An interpro description
  Arg [INTERPRO_AC]     : (optional) string An interpro accession
  Arg [TRANSLATION_ID]  : (optional) integer A translation dbID
  Arg [...]             : named arguments to FeaturePair superclass
  Example    :

    $pf =
      Bio::EnsEMBL::ProteinFeature->new( -IDESC       => $idesc,
                                         -INTERPRO_AC => $iac,
                                         @fp_args );

  Description: Instantiates a Bio::EnsEMBL::ProteinFeature
  Returntype : Bio::EnsEMBL::FeaturePair
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $proto = shift;

  my $class = ref($proto) || $proto;

  my ($idesc, $interpro_ac, $translation_id) = rearrange(['IDESC', 'INTERPRO_AC', 'TRANSLATION_ID'], @_);

  my $self = $class->SUPER::new(@_);

  # the strand of protein features is always 0
  $self->{'strand'}         = 0;
  $self->{'idesc'}          = $idesc || '';
  $self->{'interpro_ac'}    = $interpro_ac || '';
  $self->{'translation_id'} = $translation_id || '';

  return $self;
}

=head2 strand

  Arg [1]    : Ignored
  Description: Overwrites Bio::EnsEMBL::Feature->strand to not allow
             : the strand to be set.
  Returntype : int
  Status     : Stable

=cut

#do not allow the strand to be set
sub strand {
  my $self = shift;
  return $self->{'strand'};
}

=head2 idesc

  Arg [1]    : (optional) string The interpro description
  Example    : print $protein_feature->idesc();
  Description: Getter/Setter for the interpro description of this protein
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub idesc {
  my $self = shift;
  $self->{'idesc'} = shift if (@_);
  return $self->{'idesc'};
}

=head2 interpro_ac

  Arg [1]    : (optional) string The interpro accession
  Example    : print $protein_feature->interpro_ac();
  Description: Getter/Setter for the interpro accession of this protein
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub interpro_ac {
  my $self = shift;
  $self->{'interpro_ac'} = shift if (@_);
  return $self->{'interpro_ac'};
}

=head2 translation_id

  Arg [1]    : (optional) integer The dbID of the translation
  Example    : print $protein_feature->translation_id();
  Description: Getter/Setter for the translation dbID of this protein
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub translation_id {
  my $self = shift;
  $self->{'translation_id'} = shift if (@_);
  return $self->{'translation_id'};
}

1;

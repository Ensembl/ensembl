
#
# BioPerl module for ProteinFeature
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::ProteinFeature

=head1 SYNOPSIS


  my $feature = Bio::EnsEMBL::ProteinFeature->new
    (-start   => $start,
     -end     => $end,
     -hstart  => $hit_start,
     -hend    => $hit_end,
     -hname   => $hit_name);

=head1 DESCRIPTION

ProteinFeature objects represent domains or other features of interest 
on a peptide sequence.

=head1 CONTACT

  Post questions to the EnsEMBL development list ensembl-dev@sanger.ac.uk

=cut

package Bio::EnsEMBL::ProteinFeature;

use strict;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::FeaturePair);



=head2 new

  Arg [IDESC]       : (optional) string An interpro description
  Arg [INTERPRO_AC] : (optional) string An interpro accession
  Arg [...]         : named arguments to FeaturePair superclass
  Example    : $pf = Bio::EnsEMBL::ProteinFeature->new(-IDESC => $idesc,
                                                       -INTERPRO_AC => $iac,
                                                       @fp_args);
  Description: Instantiates a Bio::EnsEMBL::ProteinFeature
  Returntype : Bio::EnsEMBL::FeaturePair
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my ($idesc, $interpro_ac) = rearrange(['IDESC', 'INTERPRO_AC'], @_);

  my $self = $class->SUPER::new(@_);

  #the strand of protein features is always 0
  $self->{'strand'}      = 0;
  $self->{'idesc'}       = $idesc || '';
  $self->{'interpro_ac'} = $interpro_ac || '';

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

sub idesc{
  my $self = shift;
  $self->{'idesc'} = shift if(@_);
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

sub interpro_ac{
  my $self = shift;
  $self->{'interpro_ac'} = shift if(@_);
  return $self->{'interpro_ac'};
}


1;

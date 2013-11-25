=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::SplicingEventFeature - Object representing an alternative splicing event

=head1 SYNOPSIS

  my $ase = Bio::EnsEMBL::SplicingEventFeature->new(
    -START  => 123,
    -END    => 1045,
    -EXON_ID => $exon->dbID
  );

  # set some additional attributes
  $ase->type('flanking_exon');

=head1 DESCRIPTION

A representation of an Alternative Splicing Event Feature within the Ensembl system.

=head1 METHODS

=cut

package Bio::EnsEMBL::SplicingEventFeature;

use strict;

use POSIX;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature);



## Add to gene get_all_splicing_events



sub exon_id{
  my $self = shift;
  $self->{'exon_id'} = shift if (@_);

  if (defined $self->{'exon_id'}) {
    return $self->{'exon_id'};
  }

  return undef;
}

sub transcript_id{
  my $self = shift;
  $self->{'transcript_id'} = shift if (@_);

  if (defined $self->{'transcript_id'}) {
    return $self->{'transcript_id'};
  }

  return undef;
}

sub feature_order{
  my $self = shift;
  $self->{'feature_order'} = shift if (@_);

  if (defined $self->{'feature_order'}) {
    return $self->{'feature_order'};
  }

  return undef;
}

sub type{
  my $self = shift;
  $self->{'type'} = shift if (@_);

  if (defined $self->{'type'}) {
    return $self->{'type'};
  }

  return undef;
}

sub start{
  my $self = shift;
  $self->{'start'} = shift if (@_);

  if (defined $self->{'start'}) {
    return $self->{'start'};
  }

  return undef;
}

sub end{
  my $self = shift;
  $self->{'end'} = shift if (@_);

  if (defined $self->{'end'}) {
    return $self->{'end'};
  }

  return undef;
}


sub transcript_association{
  my $self = shift;
  $self->{'transcript_association'} = shift if (@_);

  if (defined $self->{'transcript_association'}) {
    return $self->{'transcript_association'};
  }

  return undef;
}




1;

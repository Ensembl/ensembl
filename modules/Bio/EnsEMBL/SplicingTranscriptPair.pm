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

Bio::EnsEMBL::SplicingTranscriptPair - Object representing an alternative splicing transcript pair

=head1 SYNOPSIS

  my $ase = Bio::EnsEMBL::SplicingTranscriptPair->new(
    -START  => 123,
    -END    => 1045,
    -TRANSCRIPT_ID_1 => $tran1->dbID,
    -TRANSCRIPT_ID_2 => %tran2->dbID
  );

=head1 DESCRIPTION

A representation of an Alternative Splicing Transcrript Pair within the Ensembl system.

=head1 METHODS

=cut

package Bio::EnsEMBL::SplicingTranscriptPair;

use strict;

use POSIX;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature);




sub transcript_id_1{
  my $self = shift;
  $self->{'transcript_id_1'} = shift if (@_);

  if (defined $self->{'transcript_id_1'}) {
    return $self->{'transcript_id_1'};
  }

  return undef;
}

sub transcript_id_2{
  my $self = shift;
  $self->{'transcript_id_2'} = shift if (@_);

  if (defined $self->{'transcript_id_2'}) {
    return $self->{'transcript_id_2'};
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





1;

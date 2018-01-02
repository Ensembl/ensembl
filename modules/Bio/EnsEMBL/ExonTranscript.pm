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

Bio::EnsEMBL::ExonTranscript - An Exon feature in relation to a transcript

=head1 SYNOPSIS

  use Bio::EnsEMBL::ExonTranscript;

  $feature = Bio::EnsEMBL::ExonTranscript->new(-exon => $exon, -transcript => $transcript);

=head1 DESCRIPTION

This is a feature which extends the Exon class to add
the relation to a specific transcript.

=head1 METHODS

=cut

use strict;

package Bio::EnsEMBL::ExonTranscript;

use vars qw(@ISA);

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

@ISA = qw(Bio::EnsEMBL::Exon);


=head2 new

  Arg [Exon]: The Exon object
  Arg [Transcript]: The Transcript object the exon belongs to
  Example    : $feature = Bio::EnsEMBL::ExonTranscript->new(-exon => $exon, -transcript => $transcript);
  Description: Constructs a new Bio::EnsEMBL::ExonTranscript.
  Returntype : Bio::EnsEMBL::ExonTranscript
  Exceptions : Thrown on invalid -SLICE, -STRAND arguments
  Caller     : general, subclass constructors
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  #allow this to be called as class or object method
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($exon, $transcript, $rank) = rearrange(['EXON','TRANSCRIPT', 'RANK'],@_);

  $self->{'exon'} = $exon;
  $self->{'transcript'} = $transcript;
  $self->{'rank'} = $rank;

  foreach my $attribute (keys %$exon) {
    $self->{$attribute} = $exon->{$attribute};
  }

  return $self;
}


=head2 exon

  Arg [1]    : (optional) Bio::EnsEMBL::Exon
  Example    : $exon = $exon_transcript->exon();
  Description: Getter/Setter for the exon forming this
               ExonTranscript.
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub exon {
  my $self = shift;
  $self->{'exon'} = shift if(@_);
  return $self->{'exon'};
}


=head2 transcript

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript
  Example    : print $utr->type();
  Description: Getter/Setter for the parent transcript
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transcript {
  my $self = shift;
  $self->{'transcript'} = shift if( @_ );
  return ( $self->{'transcript'} );
}

=head2 translation

    Description: Fetch the translation associated with
                 this transcript, if it exists. Return undef
                 if there is no translation, ie. a pseudogene
    Returntype : Bio::EnsEMBL::Translation or undef
    Caller     : general
    Status     : Stable

=cut

sub translation {
    my $self = shift;
    return $self->transcript()->translation();
}

=head2 get_Gene

  Description: Get the gene associated with the ExonTranscript,
               if a transcript has been set
  Returntype : Bio::EnsEMBL::Gene or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_Gene {
    my $self = shift;

    if($self->{'transcript'}) {
	return $self->{'transcript'}->get_Gene();
    }

    return;
}

=head2 rank

  Example    : print $et->rank();
  Description: Getter/Setter for the exon rank
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub rank {
  my $self = shift;
  $self->{'rank'} = shift if( @_ );
  return ( $self->{'rank'} );
}


=head2 summary_as_hash

  Example    : my $hash = $utr->summary_as_hash();
  Description: Generates a HashRef compatible with GFFSerializer. Adds
               rank and transcript attributes to the basic Feature
               hash
  Returntype : Hash
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub summary_as_hash {
  my ($self) = @_;
  my $hash = $self->SUPER::summary_as_hash();
  $hash->{'rank'} = $self->rank() if $self->rank();
  $hash->{'Parent'} = $self->transcript->display_id() if $self->transcript();
  $hash->{'source'} = $self->transcript->source() if $self->transcript();
  return $hash;
}


1;

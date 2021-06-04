=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::MANE - Object representing a MANE transcript

=head1 SYNOPSIS

  use Bio::EnsEMBL::MANE;

  $mane_transcript = Bio::EnsEMBL::MANE->new(
    -transcript    => $transcript
  );

=head1 DESCRIPTION

This is a MANE transcript
It represents an Ensembl transcript that has a matching RefSeq equivalent

=head1 METHODS

=cut

use strict;

package Bio::EnsEMBL::MANE;

use vars qw(@ISA);

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

@ISA = qw(Bio::EnsEMBL::Feature);

use constant SEQUENCE_ONTOLOGY => {
  acc  => 'SO:0000673',
  term => 'transcript',
};

=head2 new

  Arg [...]  : Named arguments passed to superclass
  Example    : $feature = Bio::EnsEMBL::MANE->new
                        (-transcript => $transcript);
  Description: Constructs a new Bio::EnsEMBL::MANE.
  Returntype : Bio::EnsEMBL::MANE
  Exceptions : 
  Caller     : general, subclass constructors
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  #allow this to be called as class or object method
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($transcript, $stable_id, $seq_region_start, $seq_region_end) = rearrange(['TRANSCRIPT', 'STABLE_ID', 'SEQ_REGION_START', 'SEQ_REGION_END'],@_);

  $self->{'transcript'} = $transcript;
  $self->{'stable_id'} = $stable_id;
  $self->{'seq_region_start'} = $seq_region_start;
  $self->{'seq_region_end'} = $seq_region_end;
  return $self;
}

=head2 transcript
  Arg [1]    : Fetch the original transcript object
  Example    : $transcript = $mane->transcript();
  Description: Getter/Setter for the transcript object
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transcript {
  my $self = shift;
  $self->{'transcript'} = shift if(@_);
  return $self->{'transcript'};
}

=head2 stable_id
  Arg [1]    : (optional) string $stable_id
  Example    : $stable_id = $mane->stable_id();
  Description: Getter/Setter for the stable_id for
               the transcript associated with this MANE transcript
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if(@_);
  return $self->{'stable_id'};
}


=head2 refseq

  Arg [1]    : Fetch the RefSeq accession associated with
               this transcript
  Example    : $refseq = $mane->refseq();
  Description: Getter/Setter for the RefSeq associated with this
               MANE transcript.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub refseq {
  my $self = shift;
  my $transcript = $self->transcript;
  my @attributes = @{ $transcript->get_all_Attributes() };
  foreach my $attribute (@attributes) {
    if ($attribute->code =~ /MANE/) {
      return $attribute->value;
    }
  }
}

=head2 type

  Arg [1]    : (optional) string mane class
  Example    : $mane = $mane->type();
  Description: Getter/Setter for the class of MANE
               associated with this transcript.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type {
  my $self = shift;
  my $transcript = $self->transcript;
  my @attributes = @{ $transcript->get_all_Attributes() };
  foreach my $attribute (@attributes) {
    if ($attribute->code =~ /MANE/) {
      return $attribute->code;
    }
  }
}

sub display_id {
  my $self = shift;
  my $id = $self->stable_id;
  return $id;
}


=head2 summary_as_hash

  Example    : my $hash = $mane->summary_as_hash();
  Description: Generates a HashRef compatible with GFFSerializer. Adds
               Gene, RefSeq accession and MANE class to the basic feature hash
  Returntype : Hash
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub summary_as_hash {
  my ($self) = @_;
  my $summary_ref = $self->SUPER::summary_as_hash;
  $summary_ref->{'refseq_match'} = $self->refseq();
  $summary_ref->{'type'} = $self->type();
  $summary_ref->{'id'} = $self->stable_id();
  my $parent = $self->transcript->get_Gene();
  $summary_ref->{'Parent'} = $parent->display_id if defined $parent;
  $summary_ref->{'strand'} = $self->transcript->strand();
  $summary_ref->{'version'} = $self->transcript->version() if $self->transcript->version();
  return $summary_ref;
}


1;

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

Bio::EnsEMBL::CDS - Object representing a CDS

=head1 SYNOPSIS

  use Bio::EnsEMBL::CDS;

  $feature = Bio::EnsEMBL::CDS->new(
    -start         => 100,
    -end           => 220,
    -strand        => -1,
    -slice         => $slice,
    -phase         => 1,
    -transcript    => $transcript
  );

=head1 DESCRIPTION

This is a CDS feature within the Ensembl CDS.
It represents the coding regions of a transcript.

=head1 METHODS

=cut

use strict;

package Bio::EnsEMBL::CDS;

use vars qw(@ISA);

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [...]  : Named arguments passed to superclass
  Example    : $feature = Bio::EnsEMBL::CDS->new
                        (-start   => 1,
                         -end     => 100,
                         -strand  => 1,
                         -slice   => $slice,
                         -dbID    => 10,
                         -transcript  => $transcript,
                         -phase   => 1);
  Description: Constructs a new Bio::EnsEMBL::CDS.
  Returntype : Bio::EnsEMBL::CDS
  Exceptions : Thrown on invalid -SLICE, -STRAND arguments
  Caller     : general, subclass constructors
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  #allow this to be called as class or object method
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($transcript, $phase, $translation_id, $seq_region_start, $seq_region_end) = rearrange(['TRANSCRIPT','PHASE', 'TRANSLATION_ID', 'SEQ_REGION_START', 'SEQ_REGION_END'],@_);

  $self->{'transcript'} = $transcript;
  $self->{'phase'} = $phase;
  $self->{'translation_id'} = $translation_id;
  $self->{'seq_region_start'} = $seq_region_start;
  $self->{'seq_region_end'} = $seq_region_end;

  return $self;
}


=head2 transcript

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript
  Example    : $transcript = $cds->transcript();
  Description: Getter/Setter for the transcript associated with this
               CDS.
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

=head2 seq_region_start

  Arg [1]    : (optional) string $seq_region_start
  Example    : $seq_region_start = $cds->seq_region_start();
  Description: Getter/Setter for the seq_region_start for this CDS.
               Overwrite default method from Feature as CDS does not have
               a table
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_start {
  my $self = shift;
  $self->{'seq_region_start'} = shift if(@_);
  return $self->{'seq_region_start'};
}

=head2 seq_region_end

  Arg [1]    : (optional) string $seq_region_end
  Example    : $seq_region_end = $cds->seq_region_end();
  Description: Getter/Setter for the seq_region_end for this CDS.
               Overwrite default method from Feature as CDS does not have
               a table
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_end {
  my $self = shift;
  $self->{'seq_region_end'} = shift if(@_);
  return $self->{'seq_region_end'};
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

=head2 translation_id

  Arg [1]    : (optional) string $translation_id
  Example    : $translation_id = $cds->translation_id();
  Description: Getter/Setter for the stable_id for the translation
               associated with this CDS.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub translation_id {
  my $self = shift;
  $self->{'translation_id'} = shift if(@_);
  return $self->{'translation_id'};
}


=head2 phase

  Arg [1]    : (optional) string $phase
  Example    : print $cds->phase();
  Description: This method returns an integer that describes
               the phase of the coding feature
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub phase {
  my $self = shift;
  $self->{'phase'} = shift if( @_ );
  return ( $self->{'phase'} || 0 );
}


=head2 summary_as_hash

  Example    : my $hash = $cds->summary_as_hash();
  Description: Generates a HashRef compatible with GFFSerializer. Adds
               Parent, source and phase to the basic feature hash
  Returntype : Hash
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub summary_as_hash {
  my ($self) = @_;
  my $hash = $self->SUPER::summary_as_hash();
  $hash->{'phase'} = $self->phase();
  my $parent = $self->transcript();
  $hash->{'Parent'} = $parent->display_id if defined $parent;
  $hash->{'source'} = $self->transcript->source() if $self->transcript();
  $hash->{'id'} = $self->translation_id() if $self->translation_id();
  $hash->{'protein_id'} = $self->translation_id() if $self->translation_id();
  return $hash;
}


1;

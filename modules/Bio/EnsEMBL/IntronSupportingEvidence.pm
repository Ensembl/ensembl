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

package Bio::EnsEMBL::IntronSupportingEvidence;

=pod


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::IntronSupportingEvidence

=head1 DESCRIPTION

Formalises an Intron with information about why it is a believed Intron. This 
serves as a parallel object to Bio::EnsEMBL::Intron which you can use 
to populate values in this field from. They are different objects though
due to Intron's non-existence as a DB data structure.

=head1 SYNOPSIS

  #Example setups a ISE from the first two Exons
  my ($five_prime_exon, $three_prime_exon) = @{$transcript->get_all_Exons()}[0..1];
  my $intron = Bio::EnsEMBL::Intron->new($five_prime_exon, $three_prime_exon);

=head1 METHODS

=cut


use strict;
use warnings;
use base qw/Bio::EnsEMBL::Feature/;

use Bio::EnsEMBL::Intron;
use Bio::EnsEMBL::Utils::Argument qw/rearrange/;
use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;

our %SUPPORTED_TYPES = map { $_ => 1 } qw/NONE DEPTH/;

=head2 new

  Arg [-ANALYSIS]     : Bio::EnsEMBL::Analysis The analysis this intron is linked to
  Arg [-START]        : int - start postion of the IntronSupportingEvidence
  Arg [-END]          : int - end position of the IntronSupportingEvidence
  Arg [-STRAND]       : int - strand the IntronSupportingEvidence is on
  Arg [-SLICE]        : Bio::EnsEMBL::Slice - the slice the IntronSupportingEvidence is on
  Arg [-INTRON]       : Bio::EnsEMBL::Intron Intron the evidence is based 
                        on. Useful if you are not specifying the location 
                        parameters as we will take them from this 
  Arg [-HIT_NAME]     : String The name of the hit
  Arg [-SCORE]        : Double The score associated with the supporting evidence
  Arg [-SCORE_TYPE]   : String The type of score we are representing
  Example             : Bio::EnsEMBL::IntronSupportingEvidence->new(
                          -ANALYSIS => $analysis, -INTRON => $intron, 
                          -SCORE => 100, -SCORE_TYPE => 'DEPTH');
  Description         : Returns a new instance of this object
  Returntype          : Bio::EnsEMBL::IntronSupportEvidence
  Exceptions          : Thrown if data is not as requested

=cut

sub new {
  my ($class, @args) = @_;
  
  my $self = $class->SUPER::new(@args);
  
  my ($intron, $hit_name, $score, $score_type, $is_splice_canonical) = 
    rearrange([qw/intron hit_name score score_type is_splice_canonical/], @args);
  
  if($intron) {
    $self->set_values_from_Intron($intron);
  }
  $self->hit_name($hit_name) if $hit_name;
  $self->score($score) if $score;
  $self->score_type($score_type) if $score_type;
  $self->is_splice_canonical($is_splice_canonical) if $is_splice_canonical;
  
  return $self;
}

=head2 set_values_from_Intron

  Arg [1]     : Bio::EnsEMBL::Intron The intron to base this object on
  Example     : $ise->set_values_from_Intron($intron);
  Description : Sets the start, end, strand and slice of this ISE instance
                using values from the given Intron object.
  Returntype  : None
  Exceptions  : Thrown if data is not as requested

=cut

sub set_values_from_Intron {
  my ($self, $intron) = @_;
  assert_ref($intron, 'Bio::EnsEMBL::Intron', 'intron');
  $self->start($intron->start());
  $self->end($intron->end());
  $self->strand($intron->strand());
  $self->slice($intron->slice());
  $self->is_splice_canonical($intron->is_splice_canonical());
  return;
}

=head2 is_splice_canonical

  Arg [1]     : Boolean
  Example     : $ise->is_splice_canonical(1);
  Description : Getter/setter for is_splice_canonical. Splice canonical 
                indicates those Introns which have a splice junction which 
                is structured as expected 
  Returntype  : Boolean
  Exceptions  : 

=cut

sub is_splice_canonical {
  my ($self, $is_splice_canonical) = @_;
  $self->{'is_splice_canonical'} = $is_splice_canonical if defined $is_splice_canonical;
  return $self->{'is_splice_canonical'};
}

=head2 get_Intron

  Arg [1]     : Bio::EnsEMBL::Transcript
  Example     : my $intron = $ise->intron($transcript);
  Description : Provides access to an Intron object by using a given transcript 
                object and its associcated array of Exons.
  Returntype  : Bio::EnsEMBL::Intron
  Exceptions  : None

=cut

sub get_Intron {
  my ($self, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');
  my $five_prime = $self->find_previous_Exon($transcript);
  my $three_prime = $self->find_next_Exon($transcript);
  return Bio::EnsEMBL::Intron->new($five_prime, $three_prime);
}

=head2 hit_name

  Arg [1]     : String name of the hit
  Example     : $ise->hit_name('hit');
  Description : Getter/setter for hit name i.e. an identifier for the alignments
  Returntype  : String
  Exceptions  : None

=cut

sub hit_name {
  my ($self, $hit_name) = @_;
  $self->{'hit_name'} = $hit_name if defined $hit_name;
  return $self->{'hit_name'};
}

=head2 score

  Arg [1]     : Number; the score associated with this feature
  Example     : $ise->score(100);
  Description : Getter/setter for score
  Returntype  : Number
  Exceptions  : None

=cut

sub score {
  my ($self, $score) = @_;
  $self->{'score'} = $score if defined $score;
  return $self->{'score'};
}

=head2 score_type

  Arg [1]     : String the enum type. Currently only allowed NONE or DEPTH
  Example     : $ise->score_type('DEPTH');
  Description : Gets and sets the type of score this instance represents
  Returntype  : String
  Exceptions  : Thrown if given an unsupported type of data

=cut

sub score_type {
  my ($self, $score_type) = @_;
  if(defined $score_type) {
    if(! $SUPPORTED_TYPES{$score_type}) {
      my $values = join(q{, }, keys %SUPPORTED_TYPES);
      throw "The score_type '$score_type' is not allowed. Allowed values are [${values}]";
    }
  }
  $self->{'score_type'} = $score_type if defined $score_type;
  return $self->{'score_type'};
}

=head2 has_linked_transcripts

  Example     : $ise->has_linked_transcripts();
  Description : Returns true if we have transcripts linked to this ISE
  Returntype  : Boolean
  Exceptions  : Thrown if we do not have an attached adaptor

=cut

sub has_linked_transcripts {
  my ($self) = @_;
  throw "No attached adaptor. Cannot find linked Transcripts unless this is a persisted object" unless $self->adaptor();
  my $transcript_ids = $self->adaptor()->list_linked_transcript_ids($self);
  return scalar(@{$transcript_ids}) ? 1 : 0;
}

=head2 equals

  Arg [1]     : Bio::EnsEMBL::IntronSupportEvidence Object to compare to
  Example     : $ise->equals($another_ise);
  Description : Asserts if the given IntronSupportEvidence instance was equal to this
  Returntype  : Boolean
  Exceptions  : None

=cut

sub equals {
  my ($self, $other) = @_;
  my $equal = $self->SUPER::equals($other);
  return 0 if ! $equal;
  return ( 
    ($self->hit_name()||q{}) eq ($other->hit_name()||q{}) &&
    ($self->score_type() eq $other->score_type()) &&
    ($self->score() == $other->score())) ? 1 : 0;
}

=head2 find_previous_Exon

  Arg [1]     : Bio::EnsEMBL::Transcript Transcript to search for the Exons from
  Example     : $ise->find_previous_Exon($transcript);
  Description : Loops through those Exons available from the Transcript and
                attempts to find one which was the 5' flanking exon. If the
                object has already been persisted we will use dbIDs to
                find the Exons
  Returntype  : Bio::EnsEMBL::Exon
  Exceptions  : None

=cut

sub find_previous_Exon {
  my ($self, $transcript) = @_;
  
  #Use DB IDs if we have them
  my $exon_id;
  if($self->adaptor()) {
    my @ids = $self->adaptor()->fetch_flanking_exon_ids($self, $transcript);
    $exon_id = $ids[0] if @ids; 
  }
  
  my $exons = $transcript->get_all_Exons();
  
  my $start = $self->start();
  my $end = $self->end();
  foreach my $exon (@{$exons}) {
    if($exon_id) {
      return $exon if $exon->dbID() == $exon_id;
      next;
    }
    if($self->strand() == 1) {
      return $exon if $exon->end() == $start-1;
    }
    else {
      return $exon if $exon->start() == $end+1;
    }
  }
  return;
}

=head2 find_next_Exon

  Arg [1]     : Bio::EnsEMBL::Transcript Transcript to search for the Exons from
  Example     : $ise->find_next_Exon($transcript);
  Description : Loops through those Exons available from the Transcript and
                attempts to find one which was the 3' flanking exon. If the
                object has already been persisted we will use dbIDs to
                find the Exons
  Returntype  : Bio::EnsEMBL::Exon
  Exceptions  : None

=cut

sub find_next_Exon {
  my ($self, $transcript) = @_;
  
  #Use DB IDs if we have them
  my $exon_id;
  if($self->adaptor()) {
    my @ids = $self->adaptor()->fetch_flanking_exon_ids($self, $transcript);
    $exon_id = $ids[1] if @ids; 
  }
  
  my $exons = $transcript->get_all_Exons();
  my $start = $self->start();
  my $end = $self->end();
  foreach my $exon (@{$exons}) {
    if($exon_id) {
      return $exon if $exon->dbID() == $exon_id;
      next;
    }
    if($self->strand() == 1) {
      return $exon if $exon->start() == $end+1;
    }
    else {
      return $exon if $exon->end() == $start-1;
    }
  }
  return;
}

1;

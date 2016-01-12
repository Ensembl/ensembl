=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Utils::IO::HtslibFileFaidx

=head1 DESCRIPTION

A class used to retrieve sequence from a FAIDX indexed FASTA file
using the htslib C bindings from Bio:DB::HTS

=cut

package Bio::EnsEMBL::Utils::IO::HtsFileFaidx;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Utils::IO::FileFaidx/;

use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::DB::HTS::Faidx;

=head2 new
  
  Arg [1]     : String; $file. Path to the FASTA file
  Arg [2]     : Boolean; $uppercase_sequence. Uppercase sequence returned from the code. Defaults to true
  Description : Builds an instance of the HtslibFaidxFasta object

=cut

sub new {
  my ($class, $file, $uppercase_sequence) = @_;
  throw 'No file given; cannot continue without one' unless $file;
  $uppercase_sequence //= 1;
  my $self = bless({}, ref($class)||$class);
  $self->file($file);
  $self->uppercase_sequence($uppercase_sequence);
  return $self;
}

sub faidx {
  my ($self) = @_;
  if(! defined $self->{faidx}) {
    $self->{faidx} = Bio::DB::HTS::Faidx->new($self->file());
  }
  return $self->{faidx};
}

=head2 uppercase_sequence
  
  Arg [1]     : Boolean; $uppercase_sequence
  Description : Controls if always uppercase sequence or not. Defaults to true

=cut

sub uppercase_sequence {
  my ($self, $uppercase_sequence) = @_;
  $self->{'uppercase_sequence'} = $uppercase_sequence if defined $uppercase_sequence;
  return $self->{'uppercase_sequence'};
}

=head2 can_access_id

  Arg[1]      : String; ID to query for
  Description : Checks the lookup to see if we have access to the id
  Returntype  : Boolean; indicating if we can query for the sequence 

=cut

sub can_access_id {
  my ($self, $id) = @_;
  my $length = $self->htslib_faidx->length($id);
  return 1 if $length > 0;
  return 0;
}

=head2 fetch_seq

  Arg[1]      : String; ID to query for
  Arg[1]      : Integer; Start of sequence to request
  Arg[1]      : Integer; Length of sequence required
  Description : Returns sequence as indexed in the underlying FAIDX store
  Returntype  : StringRef; a scalar reference to the retrieved sequence

=cut

sub fetch_seq {
  my ($self, $id, $q_start, $q_length) = @_;
  my $q_end = $q_start+$q_length;
  my $location = "$id:${q_start}-${q_end}";
  my ($seq, $length) = $self->htslib_faidx()->get_sequence($location);
  return \$seq;
}

sub DESTROY {
  my ($self) = @_;
  if(defined $self->{faidx}) {
    delete $self->{faidx};
  }
  return;
}

1;

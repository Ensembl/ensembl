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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::FastaSequenceAdaptor

=head1 DESCRIPTION

This sequence adaptor extends the BaseSequenceAdaptor code to provide an 
implementation of the .fai index lookup as defined by samtools. The code uses
this indexing system to access portions of sequence and translates Slice requests
into sensible locations for our FASTA query layer.

The adaptor must be initalised with access to a Faidx compatible object and the FASTA
file backing must use the same seq_region_name as the querying slices otherwise
we cannot return the required data.

=cut

package Bio::EnsEMBL::DBSQL::FastaSequenceAdaptor;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::DBSQL::BaseSequenceAdaptor/;

use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;
use Bio::EnsEMBL::Utils::Sequence qw/reverse_comp/;
use English qw/-no_match_vars/;
require bytes;

=head2 new
  
  Arg [1]     : FileFaidx; $faindex. A FileFaidx object or a compatible version
  Arg [2]     : Integer; $chunk_power. Size of the region to cache
  Arg [3]     : Integer; $cache_size. Number of regions to cache
  Description : Builds an instance of the FastaSequenceAdaptor

=cut

sub new {
  my ($class, $faindex, $chunk_power, $cache_size) = @_;
  my $self = $class->SUPER::new($chunk_power, $cache_size);
  $self->faindex($faindex);
  return $self;
}

=head2 fetch_by_Slice_start_end_strand
  
  Arg [1]     : Bio::EnsEMBL::Slice; $slice. Slice to fetch sequence for
  Arg [2]     : Integer; $start. Start of region to retrieve relative to the Slice (defaults to 1)
  Arg [3]     : Integer; $end. End of region to retreive relative to the Slice (defaults to length)
  Arg [4]     : Integer; $strand. Strand to fetch (defaults to 1)
  Description : Fetches sequence for the given slice. Unlike the normal SequenceAdaptor we assume
                Sequence is held in a FASTA file under the Slice's seq_region_name.
  Exception   : Thrown if we are given a circular slice
=cut

sub fetch_by_Slice_start_end_strand {
  my ( $self, $slice, $start, $end, $strand ) = @_;
  
  # Input checks
  assert_ref($slice, 'Bio::EnsEMBL::Slice', 'slice');
  if(defined $end && $start > $end && $slice->is_circular()) {
    throw "Currently we do not support circular requests";
  }
  
  #Get a new slice that spans the exact region to retrieve dna from.
  #Then constrain to seq region if it's gone negative or over the end  
  $slice = $self->expand_Slice($slice, $start, $end, $strand);
  $slice = $slice->constrain_to_seq_region();
  
  # This call is likely to barf if we try to query using a chr name
  # we do not understand. Use can_access_Slice() to make sure
  my $seq_ref = $self->_fetch_seq($slice->seq_region_name(), $slice->start(), $slice->length());
  reverse_comp($seq_ref) if $strand == -1;
  return $seq_ref;
}

=head2 faindex

  Description : Holds a reference to the Faindex object to use for sequence access

=cut

sub faindex {
  my ($self, $faindex) = @_;
  if(defined $faindex) {
    assert_ref($faindex, 'Bio::EnsEMBL::Utils::IO::FileFaidx', 'faidx');
    $self->{faindex} = $faindex;
  }
  return $self->{faindex};
}

=head2 can_access_Slice

  Description : Checks the lookup to see if we have access to the Slice given (using 
                seq region name as the ID). We reject any Circular Slice

=cut

sub can_access_Slice {
  my ($self, $slice) = @_;
  return 0 if $slice->is_circular();
  return $self->faindex()->can_access_id($slice->seq_region_name());
}

=head2 store
  
  Description : Unsupported operation. Please use a FASTA serialiser

=cut

sub store {
  throw "Unsupported operation. Cannot store sequence in a fasta file";
}

=head _fetch_raw_seq

  Description : Provides access to the underlying faindex object and returns a sequence scalar ref

=cut

sub _fetch_raw_seq {
  my ($self, $id, $start, $length) = @_;
  my $seq_ref = $self->faindex()->fetch_seq($id, $start, $length);
  return $seq_ref;
}

1;

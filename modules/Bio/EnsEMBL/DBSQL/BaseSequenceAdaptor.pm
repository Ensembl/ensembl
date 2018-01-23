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

Bio::EnsEMBL::DBSQL::BaseSequenceAdaptor

=head1 DESCRIPTION

The BaseSequenceAdaptor is responsible for the conversion of calls from 
C<fetch_by_Slice_start_end_strand()> for Sequence data into requests for a 
backing data store. In Ensembl these are the B<seqlevel> sequence region
records held in the MySQL database. 

The base adaptor also provides sequence caching based on normalisation 
technique similar to the UCSC and BAM binning indexes. The code works 
by right-shifting the requested start and end by a seq chunk power 
(by default 18 approx. 250,000bp) and then left-shifting by the same 
value. This means any value within a given window will always result in
the same value. Please see the worked examples below:

  # Equation
  p=position
  o=seq chunk power
  offset=( (p-1)>>o ) << o
  
  # Using real values
  p=1340001
  o=18
  right_shifted = (1340001-1) >> 18 == 5
  offset = 5 << 18 == 1310720

To control the size of the cache and sequences stored you can provide
the seq chunk power and the number of sequences cached.

=cut

package Bio::EnsEMBL::DBSQL::BaseSequenceAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception qw(throw);

our $SEQ_CHUNK_PWR   = 18; # 2^18 = approx. 250KB
our $SEQ_CACHE_SIZE  = 5;

=head2 fetch_by_Slice_start_end_strand

  Arg  [1]   : Bio::EnsEMBL::Slice slice
               The slice from which you want the sequence
  Arg  [2]   : Integer; $strand (optional)
               The start base pair relative to the start of the slice. Negative
               values or values greater than the length of the slice are fine.
               default = 1
  Arg  [3]   : (optional) int endBasePair
               The end base pair relative to the start of the slice. Negative
               values or values greater than the length of the slice are fine,
               but the end must be greater than or equal to the start
               count from 1
               default = the length of the slice
  Arg  [4]   : Integer; $strand (optional)
               Strand of DNA to fetch
  Returntype : StringRef (DNA requested)
  Description: Performs the fetching of DNA based upon a Slice. All fetches
               should use this method and no-other.
               
               Implementing classes are responsible for converting the
               given Slice and values into something which can be processed by 
               the underlying storage engine. Implementing class are also
               responsible for the reverse complementing of sequence.
  Exceptions : Thrown if not redefined

=cut

sub fetch_by_Slice_start_end_strand {
  my ( $self, $slice, $start, $end, $strand ) = @_;
  throw "fetch_by_Slice_start_end_strand() must be implemented";
}

=head2 can_access_Slice

  Description : Returns a boolean indiciating if the adaptor understands
                the given Slice.
  Returntype  : Boolean; if true you can get sequence for the given Slice
  Exceptions  : Thrown if not redefined

=cut

sub can_access_Slice {
  my ($self, $slice) = @_;
  throw "can_access_Slice() must be implemented";
}

=head2 expand_Slice 

  Arg  [1]    : Bio::EnsEMBL::Slice slice
                The slice from which you want the sequence
  Arg  [2]    : Integer; $strand (optional)
                The start base pair relative to the start of the slice. Negative
                values or values greater than the length of the slice are fine.
                default = 1
  Arg  [3]    : (optional) int endBasePair
                The end base pair relative to the start of the slice. Negative
                values or values greater than the length of the slice are fine,
                but the end must be greater than or equal to the start
                count from 1
                default = the length of the slice
  Arg  [4]    : Integer; $strand (optional)
                Strand of DNA to fetch
  Returntype  : Bio::EnsEMBL::Slice
  Description : Creates a new Slice which represents the requested region. Provides
                logic applicable to all SliceAdaptor instance
  Exceptions  : Thrown if the Slice is circular (we currently do not support this as generic logic)

=cut

sub expand_Slice {
  my ($self, $slice, $start, $end, $strand) = @_;
  
  $start = 1 if ! defined $start;
  $end = ($slice->end() - $slice->start()) + 1 if ! defined $end;
  $strand = 1 if ! defined $strand;

  my $right_expand = $end - $slice->length();    #negative is fine
  my $left_expand  = 1 - $start;                 #negative is fine

  my $new_slice = $slice;
  if ( $right_expand || $left_expand ) {
    $new_slice =
        $slice->strand > 0
      ? $slice->expand( $left_expand,  $right_expand )
      : $slice->expand( $right_expand, $left_expand );
  }
  
  return $new_slice;
}

=head2 new

  Arg [1]    : Int  $chunk_power; sets the size of each element of 
                    the sequence cache. Defaults to 18 which gives 
                    block sizes of ~250Kb (it is actually 2^18)
  Arg [2]    : Int  $cache_size; size of the cache. Defaults to 5 meaning
                    a cache of 1Mb if you use default values
  Example    : my $sa = $db_adaptor->get_SequenceAdaptor();
  Description: Constructor.  Calls superclass constructor and initialises
               internal cache structure.
  Returntype : Bio::EnsEMBL::DBSQL::SequenceAdaptor
  Exceptions : none
  Caller     : DBAdaptor::get_SequenceAdaptor
  Status     : Stable

=cut

sub new {
  my ($class, $chunk_power, $cache_size) = @_;
  $class = ref($class) || $class;
  my $self = bless({}, $class);
  $self->_init_seq_instance($chunk_power, $cache_size);
  return $self;
}

sub _init_seq_instance {
  my ($self, $chunk_power, $cache_size) = @_;
  $chunk_power ||= $SEQ_CHUNK_PWR;
  $cache_size ||= $SEQ_CACHE_SIZE;
  
  # use a LRU cache to limit the size
  tie my %cache, 'Bio::EnsEMBL::Utils::Cache', $cache_size, {};

  $self->{cache_size} = $cache_size;
  $self->{chunk_power} = $chunk_power;
  $self->{seq_cache_max} = ((2 ** $chunk_power) * $cache_size);
  $self->{seq_cache} = \%cache;
  return;
}

=head2 clear_cache

  Example    	: $sa->clear_cache();
  Description	: Removes all entries from the associcated sequence cache
  Returntype 	: None
  Exceptions 	: None

=cut

sub clear_cache {
  my ($self) = @_;
  %{$self->{seq_cache}} = ();
  return;
}

sub chunk_power {
  my ($self) = @_;
  return $self->{chunk_power};
}

sub cache_size {
  my ($self) = @_;
  return $self->{cache_size};
}

sub seq_cache_max {
  my ($self) = @_;
  return $self->{seq_cache_max};
}

=head2 _fetch_raw_seq

  Arg [1]     : String $id
                The identifier of the sequence to fetch.
  Arg [2]     : Integer $start
                Where to start fetching sequence from
  Arg [2]     : Integer $length
                Total length of seuqence to fetch
  Description : Performs the fetch of DNA from the backing storage 
                engine and provides it to the _fetch_seq() method
                for optional caching.
  Returntype  : ScalarRef of DNA fetched. All bases should be uppercased
  Exceptions  : Thrown if the method is not reimplemented

=cut
sub _fetch_raw_seq {
  my ($self, $id, $start, $length) = @_;
  throw "Need to implement _fetch_raw_seq(). Code must return a ScalarRef";
}

=head2 _fetch_seq

  Arg [1]     : String $id
                The identifier of the sequence to fetch.
  Arg [2]     : Integer $start
                Where to start fetching sequence from
  Arg [2]     : Integer $length
                Total length of seuqence to fetch
  Description	: If the requested region is smaller than our maximum length
                cachable region we will see if the cache already contains
                this chunk. If not we will request the region from C<_fetch_raw_seq()>
                and cache it. If the region requested is larger than 
                the maximum cacheable sequence length we pass the request
                onto C<_fetch_raw_seq()> with no caching layer.
                
                This module is also responsible for the conversion of
                requested regions into normalised region reuqests based
                on C<chunk_power>.
  Returntype 	: ScalarRef of DNA fetched. All bases should be uppercased
  Exceptions 	: Thrown when C<_fetch_raw_seq()> is not re-implemented

=cut

sub _fetch_seq {
  my ($self, $id, $start, $length) = @_;
  my $cache = $self->{'seq_cache'};
  my $cache_size = $self->cache_size();
  my $seq_chunk_pwr = $self->chunk_power();
  my $seq_cache_max = $self->seq_cache_max();

  my $seq_ref;

  # If the length of the region is less than the max size of the
  # cache then we can use the cache. The region is converted into
  # a min to max range of bins we have to query to get all the 
  # required sequence. Of these bins the cache may have all of them,
  # none of them or a combination of the two states
  if($length < $seq_cache_max) {
    my $chunk_min = ($start-1) >> $seq_chunk_pwr;
    my $chunk_max = ($start + $length - 1) >> $seq_chunk_pwr;

    # piece together sequence from cached component parts
    my $entire_seq = q{};
    for(my $i = $chunk_min; $i <= $chunk_max; $i++) {
      my $key = "${id}:${i}";
      #If it exists within the cache then add to the string. We will trim
      #down to the requested region later on
      my $cached_seq_ref = $cache->{$key};
      if($cached_seq_ref) {
        $entire_seq .= ${$cached_seq_ref};
      } 
      #Otherwise it is uncached so fetch, concat and store in the cache
      else {
        my $min = ($i << $seq_chunk_pwr) + 1;
        my $length = 1 << $seq_chunk_pwr;
        my $tmp_seq_ref = $self->_fetch_raw_seq($id, $min, $length);
        $entire_seq .= ${$tmp_seq_ref};
        $cache->{$key} = $tmp_seq_ref;
      }
    }

    # now return only the requested portion of the entire sequence
    my $min = ( $chunk_min << $seq_chunk_pwr ) + 1;
    my $seq = substr( $entire_seq, $start - $min, $length );
    $seq_ref = \$seq;
  }
  else {
    # If it was too big then just ask it from the backing store. There's
    # no use in caching this
    $seq_ref = $self->_fetch_raw_seq($id, $start, $length);
  }
  
  return $seq_ref;
}

1;

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

package XrefMapper::Methods::ChecksumBasic;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Digest::MD5;

my $DEFAULT_BATCH_SIZE = 1000;

sub new {
  my ($class, @args) = @_;
  my $self = bless({}, $class);
  
  my ($mapper, $batch_size) = rearrange([qw(mapper batch_size)], @args);
  
  throw 'No -MAPPER given' unless $mapper;
  $batch_size = $DEFAULT_BATCH_SIZE unless $batch_size;
  
  $self->mapper($mapper);
  $self->batch_size($batch_size);
  return $self;
}

sub mapper {
  my ($self, $_mapper) = @_;
  $self->{mapper} = $_mapper if defined $_mapper;
  return $self->{mapper};
}

sub batch_size {
  my ($self, $batch_size) = @_;
  $self->{batch_size} = $batch_size if defined $batch_size;
  return $self->{batch_size};
}

sub run {
  my ($self, $target, $source_id, $object_type, $db_url) = @_;
  
  my $reader = $self->_get_sequence_parser($target);
  my @results;
  my @tmp_list;
  my $batch_size = $self->batch_size();
  my $count = 0;
  while ( my $sequence = $reader->next_seq() ) {
    push(@tmp_list, $sequence);
    $count++;
    if( ($count % $batch_size) == 0) {
      my $res = $self->perform_mapping(\@tmp_list, $source_id, $object_type, $db_url);
      push(@results, @{$res});
      $self->mapper()->log_progress("Finished batch mapping of %d sequences\n", $batch_size);
      $count = 0;
      @tmp_list = ();
    }
  }
  
  #Final mapping if there were some left over
  if(@tmp_list) {
    $self->mapper()->log_progress("Finishing progess\n");
    my $res = $self->perform_mapping(\@tmp_list, $source_id, $object_type, $db_url);
    push(@results, @{$res});
    @tmp_list = ();
  }
  
  $reader->close();
  return \@results;
}

sub perform_mapping {
  my ($self, $sequences) = @_;
  throw('Override to perform the mapping you require');
}

sub _get_sequence_parser {
  my ($self, $target) = @_;
  throw "Cannot find the file '${target}'" unless -f $target;
  my $reader = Bio::SeqIO->new(-FILE => $target, -FORMAT => 'fasta');
  return $reader;
}

sub md5_checksum {
  my ($self, $sequence) = @_;
  my $digest = Digest::MD5->new();
  $digest->add($sequence->seq());
  return $digest->hexdigest();
}

1;

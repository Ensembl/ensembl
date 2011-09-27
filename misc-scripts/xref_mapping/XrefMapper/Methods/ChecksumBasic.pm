package XrefMapper::Methods::ChecksumBasic;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Digest::MD5;

my $DEFAULT_BATCH_SIZE = 10;

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
  my ($self, $target) = @_;
  
  if(! defined $target) {
    $target = $self->mapper()->core()->protein_file();
  }
  
  my $reader = $self->_get_sequence_parser($target);
  my @results;
  my @tmp_list;
  my $batch_size = $self->batch_size();
  my $count = 0;
  while ( my $sequence = $reader->next_seq() ) {
    push(@tmp_list, $sequence);
    $count++;
    if( ($count % $batch_size) == 0) {
      my $res = $self->perform_mapping(\@tmp_list);
      push(@results, @{$res});
      $self->mapper()->log_progress("Finished batch mapping of %d peptides\n", $batch_size);
      $count = 0;
      @tmp_list = ();
    }
  }
  
  #Final mapping if there were some left over
  if(@tmp_list) {
    $self->mapper()->log_progress("Finishing progess\n");
    my $res = $self->perform_mapping(\@tmp_list);
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
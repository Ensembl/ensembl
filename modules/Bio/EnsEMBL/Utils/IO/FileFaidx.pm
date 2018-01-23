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

Bio::EnsEMBL::Utils::IO::FileFaidx

=head1 DESCRIPTION

This object provides an implementation of the .fai index lookup as defined by samtools. 
This format assumes that all lines in a FASTA file are the same length (bytes and bases) 
and this information is held in a file called .fai (held in the same directory
as the FASTA file). The format is tab delimited like so:

    RecordID    RecordLength(bp)   File offset (bytes)       bp per line      bytes per line

For example:

    MT    16571   6       50      51
    1    247249719       16915   50      51

The columns are sequence ID, sequence length, offset in file, bases per line 
and bytes per line. This module will read this format (if the file is available)
or will generate it from a FASTA file. It is recommnded that you pre-generate this
by using this module with the C<write_index_to_disk()> flag on or by using the
samtools faidx binary.

Please note that we do not handle RAZF compressed indexes or FASTA files.

We recommend you use this adaptor with the fully assembled chromsome files available from 
our FTP site (primary assembly).

=cut

package Bio::EnsEMBL::Utils::IO::FileFaidx;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw/throw warning/;
use Bio::EnsEMBL::Utils::IO qw/iterate_lines work_with_file slurp/;
use English qw/-no_match_vars/;
require bytes;

=head2 new
  
  Arg [1]     : String; $file. Path to the FASTA file
  Arg [2]     : Boolean; $write_index_to_disk. Write the index to disk
  Arg [3]     : Boolean; $persist_fh. Persist the file handle between requests for sequence
  Arg [4]     : Boolean; $no_generation. Stop index generation and force the reading of an index from disk
  Arg [5]     : Boolean; $uppercase_sequence. Uppercase sequence returned from the code. Defaults to true
  Description : Builds an instance of the FaidxFasta object

=cut

sub new {
  my ($class, $fasta_file, $write_index_to_disk, $persist_fh, $no_generation, $uppercase_sequence) = @_;
  throw 'No file given; cannot continue without one' unless $fasta_file;
  $uppercase_sequence //= 1;
  my $self = bless({}, ref($class)||$class);
  $self->file($fasta_file);
  $self->write_index_to_disk($write_index_to_disk);
  $self->persist_fh($persist_fh);
  $self->no_generation($no_generation);
  $self->uppercase_sequence($uppercase_sequence);
  return $self;
}

=head2 can_access_id

  Description : Checks the lookup to see if we have access to the id

=cut

sub can_access_id {
  my ($self, $id) = @_;
  return exists $self->lookup()->{$id} ? 1 : 0;
}

=head2 file
  
  Arg [1]     : String; $file. Path to the FASTA file
  Description : Location of the FASTA file
  Exception   : Thrown if the file cannot be found

=cut

sub file {
  my ($self, $file) = @_;
  if(defined $file) {
    throw "No file found at '${file}'" unless -f $file;
    $self->{'file'} = $file;
  }
  return $self->{'file'};
}

=head2 write_index_to_disk
  
  Arg [1]     : Boolean; $write_index_to_disk. Controls if we write the index back to disk
  Description : Controls if we can write back to disk. Only run to generate the indexes

=cut

sub write_index_to_disk {
  my ($self, $write_index_to_disk) = @_;
  $self->{'write_index_to_disk'} = $write_index_to_disk if defined $write_index_to_disk;
  return $self->{'write_index_to_disk'};
}

=head2 persist_fh
  
  Arg [1]     : Boolean; $persist_fh
  Description : Controls if we leave a file handle open between sequence reads

=cut

sub persist_fh {
  my ($self, $persist_fh) = @_;
  $self->{'persist_fh'} = $persist_fh if defined $persist_fh;
  return $self->{'persist_fh'};
}

=head2 no_generation
  
  Arg [1]     : Boolean; $no_generation
  Description : Controls if we will attempt an index generation if a .fai file is missing

=cut

sub no_generation {
  my ($self, $no_generation) = @_;
  $self->{'no_generation'} = $no_generation if defined $no_generation;
  return $self->{'no_generation'};
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



=head2 index_suffix

  Description : Returns the index suffix normally used (fai)
  Returntype  : String of the index suffix

=cut

sub index_suffix {
  my ($class) = @_;
  return 'fai';
}

=head2 index_path

  Arg [1]     : String; $path. Path of the current index
  Description : Returns the index path (normally the given path plus an .fai extension)
  Returntype  : Path to the index file

=cut

sub index_path {
  my ($class, $path) = @_;
  my $suffix = $class->index_suffix();
  return "${path}.${suffix}";
}

=head2 lookup
  
  Description : Attempts to load the index from disk or create 
                it from the FASTA file in question. Once loaded it is 
                cached locally. The lookup will be written to disk
                if the write_index_to_disk attribute is true.
  Returntype  : HashRef of FASTA ID to ArrayRef of attributes
                [size, file start position, bases per line, bytes per line]
  Exception   : Thrown if we could not generate a lookup from the FASTA file or
                .fai index
=cut

sub lookup {
  my ($self) = @_;
  return $self->{lookup} if exists $self->{lookup};
  my $faindex_lookup = $self->load_from_faindex();
  if(! %{$faindex_lookup}) {
    if(!$self->no_generation()) {
      $faindex_lookup = $self->load_faindex_from_fasta();
    }
    if(! %{$faindex_lookup}) {
      throw "Cannot generate a lookup from a .fai file or from the fasta file ".$self->file();
    }
    $self->{lookup} = $faindex_lookup;
    if($self->write_index_to_disk()) {
      $self->write_faindex();
    }
  }
  else {
    $self->{lookup} = $faindex_lookup;
  }
  return $self->{lookup};
}

=head2 load_from_faindex

  Description : Loads the lookup index from a .fai file. This must be in the same location
                as the fasta file with a .fai extension. Please see samtools for more format
                information or the module description.

=cut

sub load_from_faindex {
  my ($self) = @_;
  my $index = $self->index_path($self->file());
  return {} if ! -f $index;
  my $contents = slurp($index);
  open my $fh, '<', \$contents or throw "Cannot open contents as an in-memory file: $!";
  my $lookup = $self->_load_faindex_from_fh($fh);
  close $fh;
  return $lookup;
}

=head2 _load_faindex_from_fh

  Description : Loads the .fai index from a given file handle

=cut

sub _load_faindex_from_fh {
  my ($self, $fh) = @_;
  throw "No file handle given" unless $fh;
  my %lookup;
  iterate_lines($fh, sub {
    my ($line) = @_;
    chomp $line;
    my ($id,$size,$location,$bases_per_line,$bytes_per_line) = $line =~ /^(.+) \s+ (\d+) \s+ (\d+) \s+ (\d+) \s+ (\d+) $/xms;
    # force numerification. Remember these values are from a text file
    $lookup{$id} = [$size+0, $location+0, $bases_per_line+0, $bytes_per_line+0, $id];
  });
  return \%lookup;
}

=head2 load_faindex_from_fasta

  Description : Iterates the given FASTA file looking for occurances of
                fasta record headers. We then record the start of DNA
                as a file offset, the size of the sequence, the number of
                bases per line and the number of bytes per line. This is
                identical to samtool's faidx command. Please note that
                samtools will do this a lot faster than this implementation
                as that's C and this is not.

=cut

sub load_faindex_from_fasta {
  my ($self) = @_;
  my %lookup;
  my $fasta_file = $self->file();
  my $current_values;
  work_with_file($fasta_file, 'r', sub {
    my ($fh) = @_;
    my ($line_number, $offset, $blank_line, $mismatched_lengths) = (0,0,0,0);
    while(my $line = <$fh>) {
      $line_number++;
      my $length = bytes::length($line);

      # Check for a header
      if($line =~ /^>(.+?)\s+/) {
        my $id = $1;
        # Reset the blank line and mismatched booleans since we're starting a new record
        ($blank_line, $mismatched_lengths) = (0,0);
        # Create a new current values array & add it to the hash
        # values are [sequence length, offset, bases per line, bytes per line, fasta id]
        $current_values = [0,-1,-1,-1,$id];
        $lookup{$id} = $current_values;
      }
      #Check for a blank line
      elsif(! $current_values && $line =~ /^\s+$/) {
        # Blank!
        warning "Found whitespace at line $line_number. Consider trimming";
      }
      #If not either we must be in sequence
      else {
        
        # If current record offset is set to -1 then we must set it to the current offset
        if ($current_values->[1] == -1) {
          $current_values->[1] = $offset;
        }
        
        # Minus whitespace gives us bp length
        my $bp_length = ($length-1);
        $current_values->[0] += $bp_length;
        
        # Already seen a line in the record
        if($current_values->[2] > -1) {
          # Check if we've seen a problem. Only way out of this is to continue
          # seeing blank lines until we hit another record. If not instant fail
          if($blank_line || $mismatched_lengths) {
            if($bp_length == 0) {
              # set the line number of blank line for later error reporting
              $blank_line = $line_number;
            }
            else {
              my $id = $current_values->[4];
              if($blank_line) {
                throw "FASTA record $id is misformatted. Line $blank_line is blank and embedded within a record. Please fix before rerunning this command";
              }
              my $recorded_length = $current_values->[2];
              throw "FASTA record $id is misformatted. Line $mismatched_lengths is different to the detected record length $recorded_length. Please fix before rerunning";
            }
          }
          
          # Mismatched length detection
          if($current_values->[2] != $bp_length) {
            $mismatched_lengths = $line_number;
            if($bp_length == 0) {
              $blank_line = $line_number;
            }
          }
        }
        # First line in the FASTA record
        else {
          $current_values->[2] = $bp_length; #bases per line
          $current_values->[3] = $length; #bytes per line
        }
      }
      $offset += $length;
    }
  });
  return \%lookup;
}

=head2 write_faindex

  Description : Writes the .fai index out to disk. Each line is a tab separated record recording 
                the id, size, location (file offset), bases per line and bytes per line of each
                FASTA record. This is compatible with samtools faidx. Values are stored according
                to their position in the file.

=cut

sub write_faindex {
  my ($self) = @_;
  my $lookup = $self->lookup();
  my $fasta_file = $self->file();
  my $index = $fasta_file.'.'.$self->index_suffix();
  work_with_file($index, 'w', sub {
    my ($fh) = @_;
    my @entries = sort { $a->[1] <=> $b->[1] } values %{$lookup};
    foreach my $entry (@entries) {
      my ($size, $location, $bases_per_line, $bytes_per_line, $id) = @{$entry};
      print $fh sprintf("%s\t%d\t%d\t%d\t%d\n", $id, $size, $location, $bases_per_line, $bytes_per_line);
    }
    return;
  });
  return;
}

=head2 fetch_seq

  Arg [1]     : String; $id. Identifier of the sequence in the FASTA file
  Arg [2]     : Integer; $q_start. The start of the region to find
  Arg [3]     : Integer; $q_length. The length of the region to fetch
  Description : The guts. We convert the requested start and length for an ID into
                a file position start and end. We seek to the start and read the
                requested sequence as a single operation. Line terminators are then
                substituted out and a Scalar reference handed back (de-reference to
                get to the DNA).
                
                All sequence is uppercased before returning.

=cut

sub fetch_seq {
  my ($self, $id, $q_start, $q_length) = @_;

  my $lookup = $self->lookup();
  my $info = $lookup->{$id};
  if(! $info) {
    throw "Cannot convert the $id into a valid lookup. Abort!";
  }
  my ($size, $location, $bases_per_line, $bytes_per_line) = @{$info};

  my $q_end = ($q_start + $q_length) - 1;
  # We can never request a region larger than the sequence length
  if($q_end > $size) {
    $q_end = $size;
  }
  
  my $file = $self->file();

  my $line_start = int(($q_start - 1) / $bases_per_line);
  my $line_start_position = ($q_start-1) % $bases_per_line;

  my $line_end = int(($q_end - 1) / $bases_per_line);
  my $line_end_position = ($q_end-1) % $bases_per_line;
  
  my $offset = $location + ($line_start * $bytes_per_line) + $line_start_position;
  my $end = $location + ($line_end * $bytes_per_line) + $line_end_position;
  my $length = ($end - $offset)+1;

  #Get sequence. ATMO this is a FH but why not a HTTP server in the future?
  my $seq_ref = $self->_read_from_source($file, $offset, $length);
  #Cleanup
  chomp ${$seq_ref};
  ${$seq_ref} =~ s/$INPUT_RECORD_SEPARATOR//g;
  ${$seq_ref} = uc(${$seq_ref}) if $self->uppercase_sequence();
  return $seq_ref;
}

# Open file, seek, read length and close filehandle.
# Override to read from alternative sources of FASTA formatted data
# indexed using FAIDX from sources like HTTP
sub _read_from_source {
  my ($self, $location, $offset, $length) = @_;
  my $persist_fh = $self->persist_fh();
  my $fh;
  if($persist_fh && exists $self->{fh}) {
    $fh = $self->{fh};
  }
  else {
    open $fh, '<', $location or throw "Cannot open $location for reading: $!";
    $self->{fh} = $fh if $persist_fh;
  }
  
  my $seq;
  seek($fh, $offset, 0);
  read($fh, $seq, $length);
  close $fh if ! $persist_fh;
  
  return \$seq;
}

sub DESTROY {
  my ($self) = @_;
  close $self->{fh} if $self->{fh};
}

1;

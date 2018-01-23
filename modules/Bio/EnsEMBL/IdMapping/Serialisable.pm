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

Bio::EnsEMBL::IdMapping::Serialisable - base class for serialisable objects

=head1 SYNOPSIS

  # instantiate an object which extends Serialisable
  my $object = YourObject->new(
    -DUMP_PATH  => '/tmp',
    -CACHE_FILE => 'object_cache.ser',
  );

  # serialise object to file
  my $filesize = $object->write_to_file;
  print LOG "Serialised object to file of size $filesize.\n";

  # later, create another object defining the same serialisation
  # location. specifying -LOAD_AUTO will automatically load it from the
  # serialisation file.
  my $object1 = YourObject->new(
    -DUMP_PATH  => '/tmp',
    -CACHE_FILE => 'object_cache.ser',
    -LOAD_AUTO  => 1,
  );

  # alternatively, manually load the object from file
  $object1->load_from_file;

=head1 DESCRIPTION

This is the base class for serialisable objects used by the
stable Id mapping.  It's essentially an OO wrapper for Storable,
providing a method to store (write_to_file(()) and one to retrieve
(read_from_file()) serialised objects.

This class is not instantiated itself, but rather extended by
implementing classes.

=head1 METHODS

  new
  write_to_file
  read_from_file
  dump_path
  cache_file_name
  cache_file
  loaded

=cut

package Bio::EnsEMBL::IdMapping::Serialisable;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(parse_bytes);
use Storable qw(nstore retrieve);


=head2 new

  Arg [DUMP_PATH] : String - path for object serialisation
  Arg [CACHE_FILE] : String - filename of serialised object
  Arg [AUTO_LOAD] : Boolean - determines whether object should be automatically
                    loaded on instantiation
  Description : Constructor.
  Return type : Bio::EnsEMBL::IdMapping::Serialisable implementing object
  Exceptions  : thrown on missing argument
  Caller      : implementing subclass
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dump_path, $cache_file, $auto_load) =
    rearrange([qw(DUMP_PATH CACHE_FILE AUTO_LOAD)], @_);

  throw("You must provide a cache file name") unless ($cache_file);
  
  my $self = {};
  bless ($self, $class);

  # initialise internal datastructure
  $self->{'dump_path'} = $dump_path || '.';
  $self->{'cache_file_name'} = $cache_file;

  # automatically load serialised object from file if requested
  if ($auto_load) {
    if (-s $self->cache_file) {
      $self->read_from_file;
      $self->{'loaded'} = 1;
    }
  }

  return $self;
}


=head2 write_to_file

  Example     : my $filesize = $object->write_to_file;
  Description : Serialises an object to a file (determined by
                $self->cache_file).
  Return type : String - size of serialisation file
  Exceptions  : thrown on I/O errors
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub write_to_file {
  my $self = shift;

  # create dump directory if it doesn't exist
  if (my $dump_path = $self->dump_path) {
    unless (-d $dump_path) {
      system("mkdir -p $dump_path") == 0 or
        throw("Unable to create directory $dump_path.\n");
    }
  }
  
  my $cache_file = $self->cache_file;

  eval { nstore($self->{'cache'}, $cache_file) };
  if ($@) {
    throw("Unable to store $cache_file: $@\n");
  }

  my $size = -s $cache_file;
  return parse_bytes($size);
}


=head2 read_from_file

  Example     : $object->read_from_file;
  Description : Reads a serialised object from file (determined by
                $self->cache_file).
  Return type : Bio::EnsEMBL::IdMapping::Serialisable implementing object
  Exceptions  : thrown on I/O errors
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub read_from_file {
  my $self = shift;

  my $cache_file = $self->cache_file;

  unless (-s $cache_file) {
    throw("No valid cache file found at $cache_file.");
  }

  eval { $self->{'cache'} = retrieve($cache_file); };
  if ($@) {
    throw("Unable to retrieve cache: $@");
  }

  return $self;
}


=head2 dump_path

  Arg[1]      : String - dump path for serialisation
  Example     : $object->dump_path('/tmp');
  Description : Getter/setter for the dump path for serialisation.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub dump_path {
  my $self = shift;
  $self->{'dump_path'} = shift if (@_);
  return $self->{'dump_path'};
}


=head2 cache_file_name

  Arg[1]      : String - file name for serialisation
  Example     : $object->cache_file_name('object_cache.ser');
  Description : Getter/setter for the file name for serialisation.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub cache_file_name {
  my $self = shift;
  $self->{'cache_file_name'} = shift if (@_);
  return $self->{'cache_file_name'};
}


=head2 cache_file

  Example     : my $cache_file = $object->cache_file;
  Description : Returns the path and name of the serialised object file.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub cache_file {
  my $self = shift;
  return $self->dump_path.'/'.$self->cache_file_name;
}


=head2 loaded

  Arg[1]      : Boolean - "loaded" status
  Example     : if ($object->loaded) {
                  # do something with the object that was loaded from a file
                } else {
                  # the object wasn't loaded but is new, so fill it
                }
  Description : Indicates whether a given object was loaded from its serialised
                state on disk.
  Return type : Boolean - TRUE if loaded from disk, FALSE otherwise
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub loaded {
  my $self = shift;
  $self->{'loaded'} = shift if (@_);
  return $self->{'loaded'};
}


1;


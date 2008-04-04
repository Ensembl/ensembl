package Bio::EnsEMBL::IdMapping::Serialisable;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(parse_bytes);
use Storable qw(nstore retrieve);


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


#
# getter/setters
#

sub dump_path {
  my $self = shift;
  $self->{'dump_path'} = shift if (@_);
  return $self->{'dump_path'};
}


sub cache_file_name {
  my $self = shift;
  $self->{'cache_file_name'} = shift if (@_);
  return $self->{'cache_file_name'};
}


sub cache_file {
  my $self = shift;
  return $self->dump_path.'/'.$self->cache_file_name;
}


sub loaded {
  my $self = shift;
  $self->{'loaded'} = shift if (@_);
  return $self->{'loaded'};
}


1;


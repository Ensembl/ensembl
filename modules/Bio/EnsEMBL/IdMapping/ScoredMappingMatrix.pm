package Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;

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
use Bio::EnsEMBL::IdMapping::Entry;
use Storable qw(nfreeze thaw nstore retrieve);


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dump_path, $cache_file) = rearrange([qw(DUMP_PATH CACHE_FILE)], @_);

  throw("You must provide a cache file name") unless ($cache_file);
  
  my $self = {};
  bless ($self, $class);

  # initialise internal datastructure
  $self->{'cache'}->{'matrix'} = {};
  $self->{'cache'}->{'source_list'} = {};
  $self->{'cache'}->{'target_list'} = {};
  $self->{'dump_path'} = $dump_path || '.';
  $self->{'cache_file'} = $self->dump_path."/$cache_file";

  return $self;
}


sub add_Entry {
  my $self = shift;
  my $entry = shift;

  unless ($entry and $entry->isa('Bio::EnsEMBL::IdMapping::Entry')) {
    throw("Need a Bio::EnsEMBL::IdMapping::Entry");
  }

  return $self->add_score($entry->source, $entry->target, $entry->score);
}


sub remove_Entry {
}


sub add_score {
  my $self = shift;
  my $source = shift;
  my $target = shift;
  my $score = shift;
  
  # make sure you don't put duplicates on the source and target lists
  unless (exists($self->{'cache'}->{'matrix'}->{"$source:$target"})) {
    push @{ $self->{'cache'}->{'source_list'}->{$source} }, $target;
    push @{ $self->{'cache'}->{'target_list'}->{$target} }, $source;
  }

  $self->{'cache'}->{'matrix'}->{"$source:$target"} = $score;
}


sub get_Entry {
  my $self = shift;
  my $source = shift;
  my $target = shift;

  if (exists($self->{'cache'}->{'matrix'}->{"$source:$target"})) {
    return Bio::EnsEMBL::IdMapping::Entry->new_fast(
      [$source, $target, $self->{'cache'}->{'matrix'}->{"$source:$target"}]
    );
  } else {
    return undef;
  }
}


sub get_score {
  my $self = shift;
  my $source = shift;
  my $target = shift;


  if (exists($self->{'cache'}->{'matrix'}->{"$source:$target"})) {
    return $self->{'cache'}->{'matrix'}->{"$source:$target"};
  } else {
    return undef;
  }
}


sub get_targets_for_source {
  my $self = shift;
  my $source = shift;

  return($self->{'cache'}->{'source_list'}->{$source} || []);
}


sub get_sources_for_target {
  my $self = shift;
  my $target = shift;

  return($self->{'cache'}->{'target_list'}->{$target} || []);
}


sub get_all_Entries {
  my $self = shift;

  my @result = ();
  
  foreach my $key (keys %{ $self->{'cache'}->{'matrix'} }) {
    my ($source, $target) = split(/:/, $key);
    push @result, Bio::EnsEMBL::IdMapping::Entry->new_fast(
      [$source, $target, $self->{'cache'}->{'matrix'}->{$key}]
    );
  }

  return \@result;
}


sub get_all_sources {
  return [keys %{ $_->{'source_list'} }];
}


sub get_all_targets {
  return [keys %{ $_->{'target_list'} }];
}


sub get_entry_count {
  return scalar(keys %{ $_->{'matrix'} });
}


sub get_source_count {
  return scalar(keys %{ $_->{'source_list'} });
}


sub get_target_count {
  return scalar(keys %{ $_->{'target_list'} });
}


sub get_min_max_scores {
  my $self = shift;

  my @keys = keys %{ $self->{'cache'}->{'matrix'} };

  return [undef, undef] unless (@keys);

  # initialise; this should make loop quicker
  my $min = $self->{'cache'}->{'matrix'}->{$keys[0]};
  my $max = $self->{'cache'}->{'matrix'}->{$keys[0]};
  
  foreach my $key (@keys) {
    $min = $self->{'cache'}->{'matrix'}->{$key} if ($min > $self->{'cache'}->{'matrix'}->{$key});
    $max = $self->{'cache'}->{'matrix'}->{$key} if ($max < $self->{'cache'}->{'matrix'}->{$key});
  }

  return [$min, $max];
}


sub get_average_score {
  my $self = shift;

  my @keys = keys %{ $self->{'cache'}->{'matrix'} };

  return undef unless (@keys);

  my $total = 0;
  
  foreach my $key (@keys) {
    $total += $self->{'cache'}->{'matrix'}->{$key};
  }

  return $total/scalar(@keys);
}


sub merge {
  my $self = shift;
  my $matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('You must provide a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix');
  }

  my $c = 0;

  foreach my $key (keys %{ $matrix->{'matrix'} }) {
    if (!defined($self->{'cache'}->{'matrix'}->{$key}) or
        $self->{'cache'}->{'matrix'}->{$key} < $matrix->{'matrix'}->{$key}) {
      $self->{'cache'}->{'matrix'}->{$key} = $matrix->{'matrix'}->{$key};
      $c++;
    }
  }

  return $c;
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

sub cache_file {
  my $self = shift;
  $self->{'cache_file'} = shift if (@_);
  return $self->{'cache_file'};
}


sub dump_path {
  my $self = shift;
  $self->{'dump_path'} = shift if (@_);
  return $self->{'dump_path'};
}


1;


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

use Bio::EnsEMBL::IdMapping::Serialisable;
our @ISA = qw(Bio::EnsEMBL::IdMapping::Serialisable);

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::Entry;


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  # initialise internal datastructure
  unless ($self->loaded) {
    $self->{'cache'}->{'matrix'} = {};
    $self->{'cache'}->{'source_list'} = {};
    $self->{'cache'}->{'target_list'} = {};
  }

  return $self;
}


sub flush {
  my $self = shift;

  # reset caches
  $self->{'cache'}->{'matrix'} = {};
  $self->{'cache'}->{'source_list'} = {};
  $self->{'cache'}->{'target_list'} = {};
}


sub sub_matrix {
  my $self = shift;
  my $start = shift;
  my $end = shift;

  # default to returning the full matrix if no start/end provided
  $start ||= 1;
  $end ||= $self->size;

  my $sub_matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $self->dump_path,
    -CACHE_FILE  => $self->cache_file_name,
  );
  my $i = 0;

  foreach my $key (sort keys %{ $self->{'cache'}->{'matrix'} }) {
    $i++;
    next if ($i < $start);
    last if ($i > $end);

    my ($source, $target) = split(/:/, $key);
    $sub_matrix->add_score($source, $target,
      $self->{'cache'}->{'matrix'}->{$key});
  }

  return $sub_matrix;
}


sub add_Entry {
  my $self = shift;
  my $entry = shift;

  unless ($entry and $entry->isa('Bio::EnsEMBL::IdMapping::Entry')) {
    throw("Need a Bio::EnsEMBL::IdMapping::Entry");
  }

  return $self->add_score($entry->source, $entry->target, $entry->score);
}


sub update_Entry {
  return $_[0]->add_Entry($_[1]);
}


sub remove_Entry {
  warning('Method ScoredMappingMatrix->remove_Entry not implemented (yet).');
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


sub set_score {
  my $self = shift;
  my $source = shift;
  my $target = shift;
  my $score = shift;
  
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

  return $self->{'cache'}->{'source_list'}->{$source} || [];
}


sub get_Entries_for_source {
  my $self = shift;
  my $source = shift;

  return [ map { $self->get_Entry($source, $_) }
	    @{ $self->{'cache'}->{'source_list'}->{$source} || [] } ];
}


sub get_sources_for_target {
  my $self = shift;
  my $target = shift;

  return $self->{'cache'}->{'target_list'}->{$target} || [];
}


sub get_Entries_for_target {
  my $self = shift;
  my $target = shift;

  return [ map { $self->get_Entry($_, $target) }
	    @{ $self->{'cache'}->{'target_list'}->{$target} || [] } ];
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
  my $self = shift;
  return [keys %{ $self->{'cache'}->{'source_list'} }];
}


sub get_all_targets {
  my $self = shift;
  return [keys %{ $self->{'cache'}->{'target_list'} }];
}


sub get_entry_count {
  my $self = shift;
  return scalar(keys %{ $self->{'cache'}->{'matrix'} });
}


sub size {
  return $_[0]->get_entry_count;
}


sub get_source_count {
  my $self = shift;
  return scalar(keys %{ $self->{'cache'}->{'source_list'} });
}


sub get_target_count {
  my $self = shift;
  return scalar(keys %{ $self->{'cache'}->{'target_list'} });
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

  # merge the matrices
  foreach my $key (keys %{ $matrix->{'cache'}->{'matrix'} }) {
    if (!defined($self->{'cache'}->{'matrix'}->{$key}) or
        $self->{'cache'}->{'matrix'}->{$key} < $matrix->{'cache'}->{'matrix'}->{$key}) {
      $self->{'cache'}->{'matrix'}->{$key} = $matrix->{'cache'}->{'matrix'}->{$key};
      $c++;
    }
  }

  # merge sources and target lists
  foreach my $key (keys %{ $matrix->{'cache'}->{'source_list'} }) {
    if (defined($self->{'cache'}->{'source_list'}->{$key})) {
      # need to merge lists
      my %unique = map { $_ => 1 } @{ $self->get_targets_for_source($key) };
      map { $unique{$_} = 1 } @{ $matrix->get_targets_for_source($key) };
      $self->{'cache'}->{'source_list'}->{$key} = [keys %unique];
    } else {
      # no merging needed
      $self->{'cache'}->{'source_list'}->{$key} = $matrix->{'cache'}->{'source_list'}->{$key};
    }
  }

  foreach my $key (keys %{ $matrix->{'cache'}->{'target_list'} }) {
    if (defined($self->{'cache'}->{'target_list'}->{$key})) {
      # need to merge lists
      my %unique = map { $_ => 1 } @{ $self->get_sources_for_target($key) };
      map { $unique{$_} = 1 } @{ $matrix->get_sources_for_target($key) };
      $self->{'cache'}->{'target_list'}->{$key} = [keys %unique];
    } else {
      # no merging needed
      $self->{'cache'}->{'target_list'}->{$key} = $matrix->{'cache'}->{'target_list'}->{$key};
    }
  }

  return $c;
}


sub log {
  my $self = shift;
  my $type = shift;
  my $dump_path = shift;
  
  my $debug_path = path_append($dump_path, 'debug');
  my $logfile = "$debug_path/${type}_scores.txt";
  
  open(my $fh, '>', $logfile) or
    throw("Unable to open $logfile for writing: $!");

  foreach my $entry (@{ $self->get_all_Entries }) {
    print $fh ($entry->to_string."\n");
  }

  close($fh);
}


sub to_string {
  my $self = shift;
  
  my $string = '';
  
  foreach my $entry (@{ $self->get_all_Entries }) {
    $string .= $entry->to_string."\n";
  }

  return $string;
}


1;


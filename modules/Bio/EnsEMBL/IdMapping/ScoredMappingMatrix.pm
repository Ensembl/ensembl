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

Bio::EnsEMBL::IdMapping::ScoredMappingMatrix - object holding a list of scored
Entries

=head1 SYNOPSIS

  # create a new ScoredMappingMatrix
  my $gene_scores = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH  => $dump_path,
    -CACHE_FILE => 'gene_scores.ser',
  );

  # add entries
  my $gene_scores->add_Entry($entry1);

  # serialise to file
  $gene_scores->write_to_file;

  # later, read these gene_scores from file
  my $gene_scores1 = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH  => $dump_path,
    -CACHE_FILE => 'gene_gene_scores.ser',
  );
  $gene_scores1->read_from_file;

=head1 DESCRIPTION

This object represents a collection of scores between source and target
objects.  It holds a list of Bio::EnsEMBL::IdMapping::Entry objects and
has methods to retrieve indiviual or all Entries, as well as derived
data like number of unique sources or targets, or various counts and
averages.

It is the main collection for dealing with scored relationships in the
stable Id mapping application.

=head1 METHODS

  new
  flush
  sub_matrix
  add_Entry
  update_Entry
  remove_Entry
  add_score
  set_score
  get_Entry
  get_score
  get_targets_for_source
  get_Entries_for_source
  get_sources_for_target
  get_Entries_for_target
  get_all_Entries
  get_all_sources
  get_all_targets
  get_entry_count
  size
  get_source_count
  get_target_count
  get_min_max_scores
  get_average_score
  merge
  log
  to_string

=cut

package Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::Serialisable;
our @ISA = qw(Bio::EnsEMBL::IdMapping::Serialisable);

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::IdMapping::Entry;


=head2 new

  Arg[1-N]    : see superclass
  Example     : my $gene_scores = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
                  -DUMP_PATH   => $dump_path,
                  -CACHE_FILE  => 'gene_scores.ser',
                );
  Description : Constructor.
  Return type : Bio::EnsEMBL::IdMapping::ScoredMappingMatrix
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 flush

  Example     : $gene_scores->flush;
  Description : Flushes (empties) the scoring matrix.
  Return type : none
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub flush {
  my $self = shift;

  # reset caches
  $self->{'cache'}->{'matrix'} = {};
  $self->{'cache'}->{'source_list'} = {};
  $self->{'cache'}->{'target_list'} = {};
}


=head2 sub_matrix

  Arg[1]      : Int $start - start index (inclusive)
  Arg[2]      : Int $end - end index (inclusive)
  Example     : # get the first 1000 elements in the matrix
                my $sub_matrix = $gene_scores->sub_matrix(1, 1000);
  Description : Returns a sub-matrix of the ScoredMappingMatrix. The arguments
                ($start and $end) specify the position of the first and last
                element to return (inclusive, counting starts with element 1,
                not 0)
  Return type : Bio::EnsEMBL::IdMapping::ScoredMappingMatrix
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 add_Entry

  Arg[1]      : Bio::EnsEMBL::IdMapping::Entry $entry - Entry to add
  Example     : $gene_scores->add_Entry($entry);
  Description : Adds an Entry to the scoring matrix.
  Return type : Float - the Entry's score
  Exceptions  : thrown on wrong or missing argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub add_Entry {
  my $self = shift;
  my $entry = shift;

  unless ($entry and $entry->isa('Bio::EnsEMBL::IdMapping::Entry')) {
    throw("Need a Bio::EnsEMBL::IdMapping::Entry");
  }

  return $self->add_score($entry->source, $entry->target, $entry->score);
}


=head2 update_Entry

  Arg[1]      : Bio::EnsEMBL::IdMapping::Entry $entry - Entry to update
  Example     : $gene_scores->update_Entry($entry);
  Description : Updates an Entry (or rather its score) in the scoring matrix.
                Actually delegates to add_Entry(), only there as an intuitively
                named wrapper.
  Return type : Float - the Entry's score
  Exceptions  : thrown on wrong or missing argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub update_Entry {
  return $_[0]->add_Entry($_[1]);
}


#
# not needed in the current application, so not implemented
#
sub remove_Entry {
  warning('Method ScoredMappingMatrix->remove_Entry not implemented (yet).');
}


=head2 add_score

  Arg[1]      : Int $source - source object's internal Id ("dbID")
  Arg[2]      : Int $target - target object's internal Id ("dbID")
  Arg[3]      : Float $score - score for source/target pair
  Example     : $gene_scores->add_score(1234, 5678, 0.997);
  Description : Adds a score for a source/target pair to the scoring matrix.
                This is a low-level version of add_Entry().
  Return type : Float - the score
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 set_score

  Arg[1]      : Int $source - source object's internal Id ("dbID")
  Arg[2]      : Int $target - target object's internal Id ("dbID")
  Arg[3]      : Float $score - score for source/target pair
  Example     : $gene_scores->set_score(1234, 5678, 0.997);
  Description : Sets the score for a source/target pair in the scoring matrix.
                This method is similar to add_score, but assumes that the Entry
                has been added before, so won't update the sources and target
                lists.
  Return type : Float - the score
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub set_score {
  my $self = shift;
  my $source = shift;
  my $target = shift;
  my $score = shift;
  
  $self->{'cache'}->{'matrix'}->{"$source:$target"} = $score;
}


=head2 get_Entry

  Arg[1]      : Int $source - source object's internal Id ("dbID")
  Arg[2]      : Int $target - target object's internal Id ("dbID")
  Example     : my $entry = $gene_scores->get_Entry($source_gene->id,
                  $target_gene->id);
  Description : Gets an Entry from the scoring matrix for a given source and
                target object.
  Return type : Bio::EnsEMBL::IdMapping::Entry or undef
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 get_score

  Arg[1]      : Int $source - source object's internal Id ("dbID")
  Arg[2]      : Int $target - target object's internal Id ("dbID")
  Example     : my $score = $gene_scores->get_score($source_gene->id,
                  $target_gene->id);
  Description : Gets the score from the scoring matrix for a given source and
                target object.
  Return type : Float or undef
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 get_targets_for_source

  Arg[1]      : Int $source - source object's internal Id ("dbID")
  Example     : my @targets = @{ $gene_scores->get_targets_for_source(1234) };
  Description : Returns a list of all targets which have a score against a given
                source object.
  Return type : Arrayref of Int (target objects' internal Ids)
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_targets_for_source {
  my $self = shift;
  my $source = shift;

  return $self->{'cache'}->{'source_list'}->{$source} || [];
}


=head2 get_Entries_for_source

  Arg[1]      : Int $source - source object's internal Id ("dbID")
  Example     : my @entries = @{ $gene_scores->get_Entries_for_source(1234) };
  Description : Returns a list of all Entries in the scoring matrix for a given
                source object.
  Return type : Arrayref of Bio::EnsEMBL::IdMapping::Entry objects
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_Entries_for_source {
  my $self = shift;
  my $source = shift;

  return [ map { $self->get_Entry($source, $_) }
	    @{ $self->{'cache'}->{'source_list'}->{$source} || [] } ];
}


=head2 get_sources_for_target

  Arg[1]      : Int $target - target object's internal Id ("dbID")
  Example     : my @sources = @{ $gene_scores->get_sources_for_target(5678) };
  Description : Returns a list of all sources which have a score against a given
                target object.
  Return type : Arrayref of Int (source objects' internal Ids)
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_sources_for_target {
  my $self = shift;
  my $target = shift;

  return $self->{'cache'}->{'target_list'}->{$target} || [];
}


=head2 get_Entries_for_target

  Arg[1]      : Int $target - target object's internal Id ("dbID")
  Example     : my @entries = @{ $gene_scores->get_Entries_for_target(5678) };
  Description : Returns a list of all Entries in the scoring matrix for a given
                target object.
  Return type : Arrayref of Bio::EnsEMBL::IdMapping::Entry objects
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_Entries_for_target {
  my $self = shift;
  my $target = shift;

  return [ map { $self->get_Entry($_, $target) }
	    @{ $self->{'cache'}->{'target_list'}->{$target} || [] } ];
}


=head2 get_all_Entries

  Example     : foreach my $entry (@{ $gene_scores->get_all_Entries }) {
                  # do something with the entry
                }
  Description : Returns a list of all Entries in the scoring matrix.
  Return type : Arrayref of Bio::EnsEMBL::IdMapping::Entry objects
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 get_all_sources

  Example     : my @sources = @{ $gene_scores->get_all_sources };
  Description : Returns a list of all sources in the scoring matrix.
  Return type : Arrayref of Int (source objects' internal Ids)
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_sources {
  my $self = shift;
  return [keys %{ $self->{'cache'}->{'source_list'} }];
}


=head2 get_all_targets

  Example     : my @targets = @{ $gene_scores->get_all_targets };
  Description : Returns a list of all targets in the scoring matrix.
  Return type : Arrayref of Int (target objects' internal Ids)
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_targets {
  my $self = shift;
  return [keys %{ $self->{'cache'}->{'target_list'} }];
}


=head2 get_entry_count

  Example     : my $num_entries = $gene_scores->get_entry_count;
  Description : Returns the number of Entries in the scoring matrix.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_entry_count {
  my $self = shift;
  return scalar(keys %{ $self->{'cache'}->{'matrix'} });
}


=head2 size

  Example     : my $size = $gene_scores->size;
  Description : Returns the size of the scoring matrix. Same value as returned
                by get_entry_count().
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub size {
  return $_[0]->get_entry_count;
}


=head2 get_source_count

  Example     : my $num_sources = $gene_scores->get_source_count;
  Description : Returns the number of distinct sources in the scoring matrix.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_source_count {
  my $self = shift;
  return scalar(keys %{ $self->{'cache'}->{'source_list'} });
}


=head2 get_target_count

  Example     : my $num_targets = $gene_scores->get_target_count;
  Description : Returns the number of distinct targets in the scoring matrix.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_target_count {
  my $self = shift;
  return scalar(keys %{ $self->{'cache'}->{'target_list'} });
}


=head2 get_min_max_scores

  Example     : my ($min_score, $max_score) = 
                 @{ $gene_scores->get_min_max_scores };
  Description : Returns the mininum and maximum score in the scoring matrix.
  Return type : Arrayref of Float [min_score, max_score]
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 get_average_score

  Example     : my $avg_score = $gene_scores->get_average_score;
  Description : Returns the average (mean) score in the matrix.
  Return type : Float
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 merge

  Arg[1]      : Bio::EnsEMBL::IdMapping::ScoredMappingMatrix $matrix - another
                matrix to merge with
  Example     : my $update_count = $gene_scores->merge($more_gene_scores);
  Description : Merges two scoring matrices. If there's an Entry for a
                source/target pair in both matrices, the higher score will be
                retained.
  Return type : Int - number of Entries added or updated
  Exceptions  : thrown on wrong or missing argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub merge {
  my $self   = shift;
  my $matrix = shift;

  unless ($matrix
      and $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix') )
  {
    throw(
       'You must provide a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix'
    );
  }

  my $c = 0;

  # merge the matrices
  foreach my $key ( keys %{ $matrix->{'cache'}->{'matrix'} } ) {
    if ( !defined( $self->{'cache'}->{'matrix'}->{$key} )
         or ( $self->{'cache'}->{'matrix'}->{$key} <
              $matrix->{'cache'}->{'matrix'}->{$key} ) )
    {
      $self->{'cache'}->{'matrix'}->{$key} =
        $matrix->{'cache'}->{'matrix'}->{$key};
      $c++;
    }
  }

  # merge sources and target lists
  foreach my $key ( keys %{ $matrix->{'cache'}->{'source_list'} } ) {
    if ( defined( $self->{'cache'}->{'source_list'}->{$key} ) ) {
      # need to merge lists
      my %unique =
        map { $_ => 1 } @{ $self->get_targets_for_source($key) };
      map { $unique{$_} = 1 }
        @{ $matrix->get_targets_for_source($key) };
      $self->{'cache'}->{'source_list'}->{$key} = [ keys %unique ];
    } else {
      # no merging needed
      $self->{'cache'}->{'source_list'}->{$key} =
        $matrix->{'cache'}->{'source_list'}->{$key};
    }
  }

  foreach my $key ( keys %{ $matrix->{'cache'}->{'target_list'} } ) {
    if ( defined( $self->{'cache'}->{'target_list'}->{$key} ) ) {
      # need to merge lists
      my %unique =
        map { $_ => 1 } @{ $self->get_sources_for_target($key) };
      map { $unique{$_} = 1 }
        @{ $matrix->get_sources_for_target($key) };
      $self->{'cache'}->{'target_list'}->{$key} = [ keys %unique ];
    } else {
      # no merging needed
      $self->{'cache'}->{'target_list'}->{$key} =
        $matrix->{'cache'}->{'target_list'}->{$key};
    }
  }

  return $c;
} ## end sub merge


=head2 log

  Arg[1]      : String $type - object type (e.g. 'gene')
  Arg[2]      : String $dump_path - path for writing output
  Example     : $gene_scores->log('gene', $conf->param('basedir'));
  Description : Logs all Entries in the scoring matrix to a file. Used for
                debugging.
  Return type : none
  Exceptions  : thrown on I/0 error
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 to_string

  Example     : print LOG $gene_scores->to_string, "\n";
  Description : Returns a string representation of the scoring matrix. This is
                simply a multi-line string, where each line is a stringified
                Entry.
                Useful for debugging and logging.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub to_string {
  my $self = shift;
  
  my $string = '';
  
  foreach my $entry (@{ $self->get_all_Entries }) {
    $string .= $entry->to_string."\n";
  }

  return $string;
}


1;


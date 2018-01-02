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

Bio::EnsEMBL::StableIdHistoryTree - object representing a stable ID history tree

=head1 SYNOPSIS

  my $registry = "Bio::EnsEMBL::Registry";
  my $archiveStableIdAdaptor =
    $registry->get_adaptor( 'human', 'core', 'ArchiveStableId' );

  my $stable_id = 'ENSG00000068990';
  my $history =
    $archiveStableIdAdaptor->fetch_history_tree_by_stable_id('ENSG01');

  print "Unique stable IDs in this tree:\n";
  print join( ", ", @{ $history->get_unique_stable_ids } ), "\n";

  print "\nReleases in this tree:\n";
  print join( ", ", @{ $history->get_release_display_names } ), "\n";

  print "\nCoordinates of nodes in the tree:\n\n";
  foreach my $a ( @{ $history->get_all_ArchiveStableIds } ) {
    print "  Stable ID: " . $a->stable_id . "." . $a->version . "\n";
    print "  Release: "
      . $a->release . " ("
      . $a->assembly . ", "
      . $a->db_name . ")\n";
    print "  coords: "
      . join( ', ', @{ $history->coords_by_ArchiveStableId($a) } )
      . "\n\n";
  }

=head1 DESCRIPTION

This object represents a stable ID history tree graph.

The graph is implemented as a collection of nodes (ArchiveStableId
objects) and links (StableIdEvent objects) which have positions
on an (x,y) grid. The x axis is used for releases, the y axis for
stable_ids. The idea is to create a plot similar to this (the numbers
shown on the nodes are the stable ID versions):

  ENSG001   1-------------- 2--
                                \
  ENSG003                         1-----1
                                /
  ENSG002   1-------2----------

           38      39      40    41    42

The grid coordinates of the ArchiveStableId objects in this example
would be (note that coordinates are zero-based):

  ENSG001.1               (0, 0)
  ENSG001.2               (2, 0)
  ENSG003.1 (release 41)  (3, 1) 
  ENSG003.1 (release 42)  (4, 1) 
  ENSG002.1               (0, 2)
  ENSG002.2               (1, 2)

The tree will only contain those nodes which had a change in the stable
ID version. Therefore, in the above example, in release 39 ENSG001 was
present and had version 1 (but will not be drawn there, to unclutter the
output).

The grid positions will be calculated by the API and will try to
untangle the tree (i.e. try to avoid overlapping lines).

=head1 METHODS

  new
  add_ArchiveStableIds
  add_ArchiveStableIds_for_events
  remove_ArchiveStableId
  flush_ArchiveStableIds
  add_StableIdEvents
  remove_StableIdEvent
  flush_StableIdEvents
  get_all_ArchiveStableIds
  get_all_StableIdEvents
  get_latest_StableIdEvent
  get_release_display_names
  get_release_db_names
  get_unique_stable_ids
  optimise_tree
  coords_by_ArchiveStableId
  calculate_coords
  consolidate_tree
  reset_tree
  current_dbname
  current_release
  current_assembly
  is_incomplete

=head1 RELATED MODULES

  Bio::EnsEMBL::ArchiveStableId
  Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor
  Bio::EnsEMBL::StableIdEvent

=cut

package Bio::EnsEMBL::StableIdHistoryTree;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Time::HiRes;


=head2 new

  Arg [CURRENT_DBNAME]   : (optional) name of current db
  Arg [CURRENT_RELEASE]  : (optional) current release number
  Arg [CURRENT_ASSEMBLY] : (optional) current assembly name
  Example     : my $history = Bio::EnsEMBL::StableIdHistoryTree->new;
  Description : object constructor
  Return type : Bio::EnsEMBL::StableIdHistoryTree
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = {};
  bless $self, $class;

  my ($current_dbname, $current_release, $current_assembly) =
    rearrange([qw( CURRENT_DBNAME CURRENT_RELEASE CURRENT_ASSEMBLY )], @_ );

  # initialise
  $self->{'current_dbname'} = $current_dbname;
  $self->{'current_release'} = $current_release;
  $self->{'current_assembly'} = $current_assembly;
  
  return $self;
}


=head2 add_ArchiveStableIds

  Arg[1..n]   : Bio::EnsEMBL::ArchiveStableId's @archive_ids
                The ArchiveStableIds to add to the the history tree
  Example     : my $archive_id = $archiveStableIdAdaptor->fetch_by_stable_id(
                  'ENSG00024808');
                $history->add_ArchiveStableId($archive_id);
  Description : Adds ArchiveStableIds (nodes) to the history tree. No
                calculation of grid coordinates is done at this point, you need
                to initiate this manually with calculate_coords().
                ArchiveStableIds are only added once for each release (to avoid
                duplicates).
  Return type : none
  Exceptions  : thrown on invalid or missing argument
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor::fetch_history_by_stable_id, general
  Status      : At Risk
              : under development

=cut

sub add_ArchiveStableIds {
  my ($self, @archive_ids) = @_;

  throw("You must provide one or more Bio::EnsEMBL::ArchiveStableIds to add.")
    unless (@archive_ids);

  foreach my $archive_id (@archive_ids) {
    throw("Bio::EnsEMBL::ArchiveStableId object expected.")
      unless (ref($archive_id) &&
              $archive_id->isa('Bio::EnsEMBL::ArchiveStableId'));

    $self->{'nodes'}->{$self->_node_id($archive_id)} = $archive_id;
  }
}


=head2 add_ArchiveStableIds_for_events 

  Example     : my $history = Bio::EnsEMBL::StableIdHistoryTree->new;
                $history->add_StableIdEvents($event1, $event2);
                $history->add_ArchiveStableIds_for_events;
  Description : Convenience method that adds all ArchiveStableIds for all
                StableIdEvents attached to this object to the tree.
  Return type : none
  Exceptions  : none
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor::fetch_history_by_stable_id, general
  Status      : At Risk
              : under development

=cut

sub add_ArchiveStableIds_for_events {
  my $self = shift;

  foreach my $event (@{ $self->get_all_StableIdEvents }) {
    if ($event->old_ArchiveStableId) {
      $self->add_ArchiveStableIds($event->old_ArchiveStableId);
    }
    if ($event->new_ArchiveStableId) {
      $self->add_ArchiveStableIds($event->new_ArchiveStableId);
    }
  }
}


=head2 remove_ArchiveStableId

  Arg[1]      : Bio::EnsEMBL::ArchiveStableId $archive_id
                the ArchiveStableId to remove from the tree
  Example     : $history->remove_ArchiveStableId($archive_id);
  Description : Removes an ArchiveStableId from the tree.
  Return type : none
  Exceptions  : thrown on missing or invalid argument
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor::fetch_history_by_stable_id, general
  Status      : At Risk
              : under development

=cut

sub remove_ArchiveStableId {
  my ($self, $archive_id) = @_;
    
  throw("Bio::EnsEMBL::ArchiveStableId object expected.")
    unless ($archive_id && ref($archive_id) &&
            $archive_id->isa('Bio::EnsEMBL::ArchiveStableId'));

  my %nodes = %{ $self->{'nodes'} };
  delete $nodes{$self->_node_id($archive_id)};
  $self->{'nodes'} = \%nodes;
}


=head2 flush_ArchiveStableIds

  Example     : $history->flush_ArchiveStableIds;
  Description : Remove all ArchiveStableIds from the tree.
  Return type : none
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub flush_ArchiveStableIds {
  my $self = shift;
  $self->{'nodes'} = undef;
}


#
# generate a unique node identifier
# 
sub _node_id {
  my ($self, $archive_id) = @_;
  return $archive_id->stable_id . ':' . $archive_id->db_name;
}


=head2 add_StableIdEvents 

  Arg[1..n]   : Bio::EnsEMBL::StableIdEvent's @events
                The StableIdEvents to add to the the history tree
  Example     : $history->add_StableIdEvents($event);
  Description : Adds StableIdEvents (links) to the history tree. Note that 
                ArchiveStableIds attached to the StableIdEvent aren't added to
                the tree automatically, you'll need to call
                add_ArchiveStableIds_for_events later.
  Return type : none
  Exceptions  : thrown on invalid or missing argument
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor::fetch_history_by_stable_id, general
  Status      : At Risk
              : under development

=cut

sub add_StableIdEvents {
  my ($self, @events) = @_;

  throw("You must provide one or more Bio::EnsEMBL::StableIdsEvents to add.")
    unless (@events);

  foreach my $event (@events) {
    throw("Bio::EnsEMBL::StableIdEvent object expected.")
      unless ($event->isa('Bio::EnsEMBL::StableIdEvent'));

    $self->{'links'}->{$self->_link_id($event)} = $event;
  }
}


=head2 remove_StableIdEvent 

  Arg[1]      : Bio::EnsEMBL::StableIdEvent $event
                the StableIdEvent to remove from the tree
  Example     : $history->remove_StableIdEvent($event);
  Description : Removes a StableIdEvent from the tree.
  Return type : none
  Exceptions  : thrown on missing or invalid arguments
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor::fetch_history_by_stable_id, general
  Status      : At Risk
              : under development

=cut

sub remove_StableIdEvent {
  my ($self, $event) = @_;
    
  throw("Bio::EnsEMBL::StableIdEvent object expected.") unless
    ($event && ref($event) && $event->isa('Bio::EnsEMBL::StableIdEvent'));

  my %links = %{ $self->{'links'} };
  delete $links{$self->_link_id($event)};
  $self->{'links'} = \%links;
}


=head2 flush_StableIdEvents 

  Example     : $history->flush_StableIdEvents; 
  Description : Removes all StableIdEvents from the tree.
  Return type : none
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub flush_StableIdEvents {
  my $self = shift;
  $self->{'links'} = undef; 
}


#
# generate a unique link identifier
# 
sub _link_id {
  my ($self, $event) = @_;

  my ($old_id, $old_db_name, $new_id, $new_db_name);
  if ($event->old_ArchiveStableId) {
    $old_id = $event->old_ArchiveStableId->stable_id;
    $old_db_name = $event->old_ArchiveStableId->db_name;
  }
  if ($event->new_ArchiveStableId) {
    $new_id = $event->new_ArchiveStableId->stable_id;
    $new_db_name = $event->new_ArchiveStableId->db_name;
  }

  return join(':', $old_id, $old_db_name, $new_id, $new_db_name);
}


=head2 get_all_ArchiveStableIds 

  Example     : foreach my $arch_id (@{ $history->get_all_ArchiveStableIds }) {
                  print $arch_id->stable_id, '.', $arch_id->version, "\n";
                }
  Description : Gets all ArchiveStableIds (nodes) in this tree.
  Return type : Arrayref of Bio::EnsEMBL::ArchiveStableId objects
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_ArchiveStableIds {
  my $self = shift;
  return [ values %{ $self->{'nodes'} } ]; 
}


=head2 get_all_current_ArchiveStableIds 

  Example     : foreach my $arch_id (@{ $history->get_all_current_ArchiveStableIds }) {
                  print $arch_id->stable_id, '.', $arch_id->version, "\n";
                }
  Description : Convenience method to get all current ArchiveStableIds in this
                tree.
                
                Note that no lazy loading of "current" status is done at that
                stage; as long as you retrieve your StableIdHistoryTree object
                from ArchiveStableIdAdaptor, you'll get the right answer. In
                other use cases, if you want to make sure you really get all
                current stable IDs, loop over the result of
                get_all_ArchiveStableIds() and call
                ArchiveStableId->current_version() on all of them.
  Return type : Arrayref of Bio::EnsEMBL::ArchiveStableId objects
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_current_ArchiveStableIds {
  my $self = shift;

  my @current = ();

  foreach my $arch_id (@{ $self->get_all_ArchiveStableIds }) {
    push @current, $arch_id if ($arch_id->is_current);
  }

  return \@current;
}


=head2 get_all_StableIdEvents 

  Example     : foreach my $event (@{ $history->get_all_StableIdsEvents }) {
                  print "Old stable ID: ", 
                    ($event->get_attribute('old', 'stable_id') or 'none'), "\n";
                  print "New stable ID: ", 
                    ($event->get_attribute('new', 'stable_id') or 'none'), "\n";
                  print "Mapping score: ", $event->score, "\n";
                }
  Description : Gets all StableIdsEvents (links) in this tree.
  Return type : Arrayref of Bio::EnsEMBL::StableIdEvent objects
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_StableIdEvents {
  my $self = shift;
  return [ values %{ $self->{'links'} } ]; 
}


=head2 get_latest_StableIdEvent

  Arg[1]      : Bio::EnsEMBL::ArchiveStableId $arch_id - the stable ID to get
                the latest Event for
  Example     : my $arch_id = Bio::EnsEMBL::ArchiveStableId->new(
                  -stable_id => 'ENSG00001'
                );
                my $event = $history->get_latest_Event($arch_id);
  Description : Returns the latest StableIdEvent found in the tree where a given
                stable ID is the new stable ID. If more than one is found (e.g.
                in a merge scenario in the latest mapping), preference is given
                to self-events.
  Return type : Bio::EnsEMBL::StableIdEvent
  Exceptions  : thrown on missing or wrong argument
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor::add_all_current_to_history, general
  Status      : At Risk
              : under development

=cut

sub get_latest_StableIdEvent {
  my $self = shift;
  my $arch_id = shift;
  
  unless ($arch_id and $arch_id->isa('Bio::EnsEMBL::ArchiveStableId')) {
    throw("Need a Bio::EnsEMBL::ArchiveStableId.");
  }

  my @all_events = @{ $self->get_all_StableIdEvents };
  my @self_events = ();

  while (my $event = shift(@all_events)) {
    if ($event->new_ArchiveStableId and
        $event->new_ArchiveStableId->stable_id eq $arch_id->stable_id) {
      push @self_events, $event;
    }
  }

  my @sorted = sort { $b->new_ArchiveStableId->release <=>
                      $a->new_ArchiveStableId->release } @self_events;
  
  # give priority to self events
  my $latest;
  while ($latest = shift @sorted) {
    last if (($latest->old_ArchiveStableId and
              $latest->old_ArchiveStableId->stable_id eq $arch_id->stable_id)
             or !$latest->old_ArchiveStableId);
  }

  return $latest;
}


=head2 get_release_display_names

  Example     : print "Unique release display_names in this tree:\n"
                foreach my $name (@{ $history->get_release_display_names }) {
                  print "  $name\n";
                }
  Description : Returns a chronologically sorted list of unique release
                display_names in this tree.

                This method can be used to determine the number of columns when
                plotting the history tree.
  Return type : Arrayref of strings.
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_release_display_names {
  my $self = shift;
  
  my @display_names = map { $_->[1] } @{ $self->_sort_releases };

  return \@display_names;
}


=head2 get_release_db_names

  Example     : print "Unique release db_names in this tree:\n"
                foreach my $name (@{ $history->get_release_db_names }) {
                  print "  $name\n";
                }
  Description : Returns a chronologically sorted list of unique release
                db_names in this tree.
  Return type : Arrayref of strings.
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_release_db_names {
  my $self = shift;
  
  my @db_names = map { $_->[0] } @{ $self->_sort_releases };

  return \@db_names;
}


#
# Create a chronologically sorted list of releases.
#
# Return type : Arrayref of arrayrefs (db_name, release)
#
sub _sort_releases {
  my $self = shift;

  unless ($self->{'sorted_tree'}->{'releases'}) {

    my %unique = ();

    foreach my $archive_id (@{ $self->get_all_ArchiveStableIds }) {
      $unique{join(':', $archive_id->db_name, $archive_id->release)} = 1;
    }

    # sort releases by release number, then db_name; this should get them into
    # chronological order
    my @releases = sort { $a->[1] <=> $b->[1] || $a->[0] cmp $b->[0] }
      map { [ split(/:/, $_) ] } keys(%unique);

    $self->{'sorted_tree'}->{'releases'} = \@releases;
  
  }

  return $self->{'sorted_tree'}->{'releases'};
}


=head2 get_unique_stable_ids 

  Example     : print "Unique stable IDs in this tree:\n"
                foreach my $id (@{ $history->get_unique_stable_ids }) {
                  print "  $id\n";
                }
  Description : Returns a list of unique stable IDs in this tree. Version is not
                taken into account here. This method can be used to determine
                the number of rows when plotting the history with each stable ID
                occupying one line.

                Sort algorithm will depend on what was chosen when the sorted
                tree was generated. This ranges from a simple alphanumeric sort
                to algorithms trying to untangle the history tree. If no
                pre-sorted data is found, an alphanumerically sorted list will
                be returned by default.
  Return type : Arrayref of strings.
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_unique_stable_ids {
  my $self = shift;
  
  unless ($self->{'sorted_tree'}->{'stable_ids'}) {
    $self->{'sorted_tree'}->{'stable_ids'} = $self->_sort_stable_ids;
  }
  
  return $self->{'sorted_tree'}->{'stable_ids'};
}


#
# Returns a list of stable IDs in this history tree, sorted alphabetically.
# This is the simplest sort function used and doesn't try to untangle the tree.
#
# Return type : Arrayref
#
sub _sort_stable_ids {
  my $self = shift;
  my %unique = map { $_->stable_id => 1 } @{ $self->get_all_ArchiveStableIds };
  return [sort keys %unique];
}


=head2 optimise_tree

  Arg [1]     : (optional) Float $time_limit
                Optimise tree normally runs until it hits a minimised state
                but this can take a very long time. Therefore you can
                opt to bail out of the optimisation early. Specify the
                time in seconds. Floating point values are supported should you
                require sub-second limits              
  Example     : $history->optimise_tree;
  Description : This method sorts the history tree so that the number of
                overlapping branches is minimised (thus "untangling" the tree).
                
                It uses a clustering algorithm for this which iteratively moves
                the nodes with the largest vertical distance next to each other
                and looking for a mininum in total branch length. This might not
                produce the overall optimum but usually converges on a local
                optimum very quickly.
  Return type : none
  Exceptions  : none
  Caller      : calculate_coords
  Status      : At Risk
              : under development

=cut

sub optimise_tree {
  my $self = shift;
  my $time_limit = shift;

  # get all non-self events
  my @links;
  foreach my $event (@{ $self->get_all_StableIdEvents }) {
    next unless ($event->old_ArchiveStableId and $event->new_ArchiveStableId);
    my $old_id = $event->old_ArchiveStableId->stable_id;
    my $new_id = $event->new_ArchiveStableId->stable_id;
    push @links, [$old_id, $new_id] if ($old_id ne $new_id);
  }

  # get initial list of sorted unique stable IDs and put them into a position
  # lookup hash
  my $i = 0;
  my %pos = map { $_ => $i++ } @{ $self->_sort_stable_ids };

  my $opt_length;
  my $successive_fails = 0;
  my $k = 0;
  my %seen;

  # for debug purposes:
  # find the number of permutations for the given number of stable IDs
  my $fact = $self->_factorial(scalar(keys %pos));

  my $starting_time = Time::HiRes::time();

  OPT:
  while ($successive_fails < 100) {

    if(defined $time_limit) {
      my $current_time = Time::HiRes::time();
      my $diff = $current_time - $starting_time;
      last OPT if $diff > $time_limit;
    }

    # sort links by vertical distance
    #warn "sorting\n";
    $self->_sort_links(\@links, \%pos);

    # loop over sorted links
    SORTED:
    foreach my $link (@links) {
      
      #warn "  trying ".join('-', @$link)."\n";

      $k++;
      
      # remember last sort order
      my %last = %pos;
      
      #my $this_order = join(':', sort { $pos{$a} <=> $pos{$b} } keys %pos);
      #warn "    before $this_order\n";

      # try both to move bottom node next to top node's current position and
      # top node next to bottom node's position - one of the methods might give
      # you better results
      DIRECT:
      foreach my $direction (qw(up down)) {

        # move the nodes next to each other
        $self->_move_nodes($link, \%pos, $direction);

        # next if we've seen this sort order before
        my $new_order = join(':', sort { $pos{$a} <=> $pos{$b} } keys %pos);
        #warn "    after ($direction) $new_order\n";
        if ($seen{$new_order}) {
          #warn "      seen\n";
          %pos = %last;
          next DIRECT;
        }
        $seen{$new_order} = 1;

        # calculate total link length for this sort order
        my $total_length = $self->_total_link_length(\@links, \%pos);

        if (!$opt_length or $total_length < $opt_length) {
          #warn "      better ($total_length/$opt_length)\n";
          $opt_length = $total_length;
          $successive_fails = 0;
          next OPT;
        } else {
          #warn "      worse ($total_length/$opt_length)\n";
          %pos = %last;
          $successive_fails++;
        }
      }
      
    }

    last OPT;
    
  }

  #warn "Needed $k tries (of $fact) to find optimal tree.\n";

  my @best = sort { $pos{$a} <=> $pos{$b} } keys %pos;
  $self->{'sorted_tree'}->{'stable_ids'} = \@best;
}


#
# find the number of permutations for a give array size.
# used for debugging code (compare implemented algorithm to looping over all
# possible permutations).
#
sub _factorial {
  my ($self, $n) = @_;
  my $s = 1;
  $s *= $n-- while $n > 0;
  return $s;
}


#
# sort links by vertical distance
#
sub _sort_links {
  my ($self, $links, $pos) = @_;

  my @lookup;

  foreach my $link (@$links) {
    my $dist = $pos->{$link->[0]} - $pos->{$link->[1]};
    $dist = -$dist if ($dist < 0);
    push @lookup, [$dist, $link];
    #warn " $dist ".join(' ', @$link)."\n";
  }

  @$links = map { $_->[1] } sort { $b->[0] <=> $a->[0] } @lookup;
}


#
# make two nodes adjacent by moving the second node next to the first node
# all other node coordinates are adjusted accordingly
#
sub _move_nodes {
  my ($self, $link, $pos, $direction) = @_;

  my $first_pos = $pos->{$link->[0]};
  my $second_pos = $pos->{$link->[1]};

  # swap positions if necessary
  if ($first_pos > $second_pos) {
    my $tmp = $second_pos;
    $second_pos = $first_pos;
    $first_pos = $tmp;
  }
  #warn "      $first_pos:$second_pos\n";

  foreach my $p (keys %$pos) {
    
    my $val = $pos->{$p};
    
    #warn "      $p $val\n";
    if ($direction eq 'up') {
      if ($val > $first_pos and $val < $second_pos) {
        $val++;
      } elsif ($val == $second_pos) {
        $val = $first_pos + 1;
      }
    } else {
      if ($val > $first_pos and $val < $second_pos) {
        $val--;
      } elsif ($val == $first_pos) {
        $val = $second_pos - 1;
      }
    }
    
    #warn "      $p $val\n";
    $pos->{$p} = $val;
    #warn "\n";
  }
}


#
# calculate the total link (vertical distance) length based on this sort order
#
sub _total_link_length {
  my ($self, $links, $pos) = @_;

  my $total_length;

  foreach my $link (@$links) {
    my $length = $pos->{$link->[0]} - $pos->{$link->[1]};
    $length = -$length if ($length < 0);
    $total_length += $length;
  }

  return $total_length;
}


=head2 coords_by_ArchiveStableId 

  Arg[1]      : Bio::EnsEMBL::ArchiveStableId $archive_id
                The ArchiveStableId to get tree grid coordinates for
  Example     : my ($x, $y) =
                  @{ $history->coords_by_ArchiveStableId($archive_id) };
                print $archive_id->stable_id, " coords: $x, $y\n";
  Description : Returns the coordinates of an ArchiveStableId in the history
                tree grid. If the ArchiveStableId isn't found in this tree, an
                empty list is returned.
                
                Coordinates are zero-based (i.e. the top leftmost element in
                the grid has coordinates [0, 0], not [1, 1]). This is to
                facilitate using them to create a matrix as a two-dimensional
                array of arrays.
  Return type : Arrayref (x coordinate, y coordinate)
  Exceptions  : thrown on wrong argument type
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub coords_by_ArchiveStableId {
  my ($self, $archive_id) = @_;

  throw("Bio::EnsEMBL::ArchiveStableId object expected.")
    unless ($archive_id and ref($archive_id) and
      $archive_id->isa('Bio::EnsEMBL::ArchiveStableId'));
  
  return $self->{'sorted_tree'}->{'coords'}->{$self->_node_id($archive_id)}
    || [];
}


=head2 calculate_coords

  Arg [1]     : (optional) Float $time_limit
                Optimise tree normally runs until it hits a minimised state
                but this can take a very long time. Therefore you can
                opt to bail out of the optimisation early. Specify the
                time in seconds. Floating point values are supported should you
                require sub-second limits
  Example     : $history->calculate_coords;
  Description : Pre-calculates the grid coordinates of all nodes in the tree.
  Return type : none
  Exceptions  : none
  Caller      : ArchiveStableIdAdaptor::fetch_history_by_stable_id
  Status      : At Risk
              : under development

=cut

sub calculate_coords {
  my $self = shift;
  my $time_limit = shift;

  # reset any previous tree cordinate calculations
  $self->reset_tree;

  # the "master" information for the sorted tree is stored as the sorted lists
  # of releases (x) and stable IDs (y). Sort them now.
  my $db_names = $self->get_release_db_names;

  # untangle tree by sorting stable IDs appropriately
  $self->optimise_tree($time_limit);
  my $stable_ids = $self->get_unique_stable_ids;
  
  # for performance reasons, additionally store coordinates in a lookup hash
  foreach my $archive_id (@{ $self->get_all_ArchiveStableIds }) {
  
    # coordinates are positions in the sorted lists
    my $x = $self->_index_of($archive_id->db_name, $db_names);
    my $y = $self->_index_of($archive_id->stable_id, $stable_ids);
  
    $self->{'sorted_tree'}->{'coords'}->{$self->_node_id($archive_id)} =
      [ $x, $y ];
  }
}

#
# Description : Returns the index of an element in an array
# Example     : my @array = (a, b, c);
#               my $i = _index_of('b', \@array); # will return 1
# Return type : Int (or undef if element is not found in array)
#
sub _index_of {
  my ($self, $element, $arrayref) = @_;

  throw("Expecting arrayref argument.") unless (ref($arrayref) eq 'ARRAY');

  my @array = @$arrayref;

  while (my $e = pop(@array)) {
    return scalar(@array) if ($e eq $element);
  }

  return undef;
}


=head2 consolidate_tree

  Example     : $history->consolidate_tree;
  Description : Consolidate the history tree. This means removing nodes where
                there wasn't a change and bridging gaps in the history. The end
                result will be a sparse tree which only contains the necessary
                information.
  Return type : none
  Exceptions  : none
  Caller      : ArchiveStableIdAdaptor->fetch_history_tree_by_stable_id
  Status      : At Risk
              : under development

=cut

sub consolidate_tree {
  my $self = shift;

  #
  # get all self-events and creations/deletions and sort them (by stable ID and
  # chronologically)
  #
  my @event_lookup;
  
  foreach my $event (@{ $self->get_all_StableIdEvents }) {

    my $old_id = $event->old_ArchiveStableId;
    my $new_id = $event->new_ArchiveStableId;

    if (!$old_id or !$new_id or ($old_id->stable_id eq $new_id->stable_id)) {
      if ($old_id) {
        push @event_lookup, [$old_id->stable_id, $old_id->release, 
          $old_id->db_name, $event];
      } else {
        push @event_lookup, [$new_id->stable_id, $new_id->release - 1,
          $new_id->db_name, $event];
      }
    }
  }

  my @self_events = map { $_->[3] }
    sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] cmp $b->[2] }
      @event_lookup;

  #
  # consolidate tree
  #
  my $last = shift(@self_events);

  while (my $event = shift(@self_events)) {

    my $lo = $last->old_ArchiveStableId;
    my $ln = $last->new_ArchiveStableId;
    my $eo = $event->old_ArchiveStableId;
    my $en = $event->new_ArchiveStableId;

    if ($lo and $eo and $en and $lo->stable_id eq $eo->stable_id
        and $lo->version eq $eo->version) {

      # this removes redundant nodes and connects the remaining nodes:
      #
      # o--o--o  ->  o-----o
      # 1  1  1      1     1

      #warn 'A: '.$last->ident_string.' | '.$event->ident_string."\n";

      $self->remove_StableIdEvent($last);
      $self->remove_StableIdEvent($event);

      $event->old_ArchiveStableId($lo);

      $self->add_StableIdEvents($event);

    } elsif ($ln and $eo and $ln->db_name ne $eo->db_name
        and $ln->stable_id eq $eo->stable_id and $ln->version eq $eo->version) {
        
      # try to brigde gaps

      if ($en) {
        
        # o--o  o--o  ->  o--o-----o
        # 1  2  2  2      1  2     2
        #
        #    o  o--o  ->  o-----o
        #    1  1  1      1     1
        
        #warn 'X: '.$last->ident_string.' | '.$event->ident_string."\n";

        $self->remove_StableIdEvent($event);
        $event->old_ArchiveStableId($ln);
        $self->add_StableIdEvents($event);

      } elsif ($lo) {
        
        # there's a deletion event, deal with it differently

        if ($lo->version eq $ln->version) {
        
          # o--o  o  ->  o-----o
          # 1  1  1      1     1
          
          #warn 'Y: '.$last->ident_string.' | '.$event->ident_string."\n";

          $self->remove_StableIdEvent($last);
          $last->new_ArchiveStableId($eo);
          $self->add_StableIdEvents($last);

        } else {

          # o--o  o  ->  o--o--o
          # 1  2  2      1  2  2
          
          #warn 'Z: '.$last->ident_string.' | '.$event->ident_string."\n";

          $self->remove_StableIdEvent($event);
          $event->old_ArchiveStableId($ln);
          $event->new_ArchiveStableId($eo);
          $self->add_StableIdEvents($event);

        }

      } else {

        # creation followed by deletion in next mapping
        #
        # o  o  ->  o--o
        # 1  1      1  1

        #warn 'Q: '.$last->ident_string.' | '.$event->ident_string."\n";

        $self->remove_StableIdEvent($last);
        $self->remove_StableIdEvent($event);
        $event->old_ArchiveStableId($ln);
        $event->new_ArchiveStableId($eo);
        $self->add_StableIdEvents($event);

      }

    } else {
      #warn 'C: '.$last->ident_string.' | '.$event->ident_string."\n";
    }
  
    $last = $event;
  }
  
  # now add ArchiveStableIds of the remaining events to the tree
  $self->add_ArchiveStableIds_for_events;
}


=head2 reset_tree

  Example     : $history->reset_tree;
  Description : Resets all pre-calculated tree grid data. Mostly used internally
                by methods that modify the tree.
  Return type : none
  Exceptions  : none
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub reset_tree {
  my $self = shift;
  $self->{'sorted_tree'} = undef;
}


=head2 current_dbname

  Arg[1]      : (optional) String $dbname - the dbname to set
  Example     : my $dbname = $history->current_dbname;
  Description : Getter/setter for current dbname.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub current_dbname {
  my $self = shift;
  $self->{'current_dbname'} = shift if (@_);
  return $self->{'current_dbname'};
}


=head2 current_release

  Arg[1]      : (optional) Int $release - the release to set
  Example     : my $release = $history->current_release;
  Description : Getter/setter for current release.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub current_release {
  my $self = shift;
  $self->{'current_release'} = shift if (@_);
  return $self->{'current_release'};
}


=head2 current_assembly

  Arg[1]      : (optional) String $assembly - the assembly to set
  Example     : my $assembly = $history->current_assembly;
  Description : Getter/setter for current assembly.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub current_assembly {
  my $self = shift;
  $self->{'current_assembly'} = shift if (@_);
  return $self->{'current_assembly'};
}


=head2 is_incomplete

  Arg[1]      : (optional) Boolean $incomplete 
  Example     : if ($history->is_incomplete) {
                  print "Returned tree is incomplete due to too many mappings
                    in the database.\n";
                }
  Description : Getter/setter for incomplete flag. This is used by
                ArchiveStableIdAdaptor to indicate that it finished building
                the tree prematurely due to too many mappins in the db and can
                be used by applications to print warning messages.
  Return type : Boolean
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_incomplete {
  my $self = shift;
  $self->{'incomplete'} = shift if (@_);
  return $self->{'incomplete'};
}


1;


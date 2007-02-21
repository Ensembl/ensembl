package Bio::EnsEMBL::StableIdHistoryTree;

=head1 NAME

Bio::EnsEMBL::StableIdHistoryTree - object representing a stable ID history tree

=head1 SYNOPSIS


=head1 DESCRIPTION

This object represents a stable ID history tree graph.

The graph is implemented as a collection of nodes (ArchiveStableId objects) and
links (StableIdEvent objects) which have positions on an (x,y) grid. The x axis
is used for releases, the y axis for stable_ids. The idea is to create a plot
similar to this (the numbers shown on the nodes are the stable ID versions):

ENSG001   1-------------- 2--
                              \
ENSG003                         1-----1
                              /
ENSG002   1-------2----------

         38      39      40    41    42

The grid coordinates of the ArchiveStableId objects in this example would
be (note that coordinates are zero-based):

ENSG001.1               (0, 0)
ENSG001.2               (2, 0)
ENSG003.1 (release 41)  (3, 1) 
ENSG003.1 (release 42)  (4, 1) 
ENSG002.1               (0, 2)
ENSG002.2               (1, 2)

The tree will only contain those nodes which contain a change in the stable
ID version. Therefore, in the above example, in release 39 ENSG001 was
present and had version 1 (but will not be drawn there, to unclutter the
output).

The grid positions will be calculated by the API and will ideally make sure
you don't get overlapping lines (not fully implemented yet).

=head1 METHODS

add_ArchiveStableIds
remove_ArchiveStableId
add_StableIdEvents
remove_StableIdEvent
get_all_ArchiveStableIds
get_all_StableIdEvents
get_release_display_names
get_release_db_names
get_unique_stable_ids
coords_by_ArchiveStableId
calculate_simple_coords
reset_tree

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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 new

  Example     : my $history_tree = Bio::EnsEMBL::StableIdHistoryTree->new;
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
  
  return $self;
}


=head2 add_ArchiveStableIds

  Arg[1]      : Bio::EnsEMBL::ArchiveStableId's @archive_ids
                The ArchiveStableIds to add to the the history tree
  Example     : my $archive_id = $archiveStableIdAdaptor->fetch_by_stable_id(
                  'ENSG00024808');
                $history_tree->add_ArchiveStableId($archive_id);
  Description : Adds ArchiveStableIds (nodes) to the history tree. No
                calculation of grid coordinates is done, you need to initiate
                this manually with calculate_coords(). ArchiveStableIds are only
                added once for each release (to avoid duplicates).
  Return type : none
  Exceptions  : thrown on invalid or missing argument
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor::fetch_history_tree_by_stable_id, general
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

  # reset pre-calculated tree data (needs to be redone after adding objects)
  $self->reset_tree;
}


=head2 remove_ArchiveStableId

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
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
  Example     : $history_tree->add_StableIdEvent($event);
  Description : Adds StableIdEvents (links) to the history tree.
  Return type : none
  Exceptions  : thrown on invalid or missing argument
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor::fetch_history_tree_by_stable_id, general
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

    # also add ArchiveStableIds linked to this event
    if ($event->old_ArchiveStableId) {
      $self->add_ArchiveStableIds($event->old_ArchiveStableId);
    }
    if ($event->new_ArchiveStableId) {
      $self->add_ArchiveStableIds($event->new_ArchiveStableId);
    }
  }

  # reset pre-calculated tree data (needs to be redone after adding objects)
  $self->reset_tree;
}


=head2 remove_StableIdEvent 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      : At Risk
              : under development

=cut

sub remove_StableIdEvent {
  my ($self, $event) = @_;
    
  throw("Bio::EnsEMBL::StableIdEvent object expected.") unless
    ($event && ref($event) && $event->isa('Bio::EnsEMBL::StableIdEvent'));

  my %links = %{ $self->{'links'} };
  delete $links{$self->_link_id($event)};
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

  Example     : foreach my $arch_id (@{ $history_tree->get_all_ArchiveStableIds }) {
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


=head2 get_all_StableIdEvents 

  Example     : foreach my $event (@{ $history_tree->get_all_StableIdsEvents }) {
                  print "Old stable ID: ", 
                    $event->old_ArchiveStableId->stable_id, "\n";
                  print "New stable ID: ", 
                    $event->new_ArchiveStableId->stable_id, "\n";
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


=head2 get_release_display_names

  Example     : print "Unique release display_names in this tree:\n"
                foreach my $name (@{ $history_tree->get_release_display_names }) {
                  print "  $name\n";
                }
  Description : Returns a chronologically sorted list of unique release
                display_names in this tree.

                Note that these display_names will differ from the value of
                ArchiveStableId->release in some rare cases where the release
                name stored in the database is not unique. In these cases, a
                "version number" is appended to the release number to create
                the display_name.
                
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
                foreach my $name (@{ $history_tree->get_release_db_names }) {
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


=head2 _sort_releases

  Example     : 
  Description : Create a chronologically sorted list of releases.

                Before release 21, sometimes several releases had the same
                number (because this number indicated schema version then which
                wasn't changed in every release). To get unique release
                identifiers we therefore need to sort this out by adding
                "version numbers" to the release.
  Return type : Arrayref of arrayrefs (db_name, release)
  Exceptions  : none
  Caller      : get_release_display_names
  Status      : At Risk
              : under development

=cut

sub _sort_releases {
  my $self = shift;

  unless ($self->{'sorted_tree'}->{'releases'}) {

    my @archive_ids = @{ $self->get_all_ArchiveStableIds };
    my $r;
    my @releases;

    # get all releases and their associated database names.
    # the combination of release number and database name is unique, so we will
    # use this to identify our releases
    while (my $archive_id = shift(@archive_ids)) {
      $r->{$archive_id->release}->{$archive_id->db_name} = 1;
    }

    foreach my $release (keys %$r) {
      
      my @db_names = sort keys %{ $r->{$release} };
      
      if (scalar(@db_names) > 1) {
        # more than one release with this number.
        # we need to create multiple display_names
        my $i = 0;
        foreach my $db_name (@db_names) {
          my $name = "$release." . ++$i;
          push @releases, [ $db_name, $name ];
        }

      } else {
        # just a single release with this number, use it directly
        push @releases, [ $db_names[0], $release ];
      }
    }

    # sort releases by release number, then db_name; this should get them into
    # chronological order
    @releases = sort { $a->[1] <=> $b->[1] || $a->[0] cmp $b->[0] } @releases;

    $self->{'sorted_tree'}->{'releases'} = \@releases;
  
  }

  return $self->{'sorted_tree'}->{'releases'};
}


=head2 get_unique_stable_ids 

  Example     : print "Unique stable IDs in this tree:\n"
                foreach my $id (@{ $history_tree->get_unique_stable_ids }) {
                  print "  $id\n";
                }
  Description : Returns a list of unique stable IDs in this tree. Version is not
                taken into account here. This method can be used to determine
                the number of rows when plotting the history with each stable ID
                occupying one line.

                Sort algorithm will depend on what was chosen when the sorted
                tree was generated. This ranges from a simple alphanumeric sort
                to algorithms trying to disentangle the history tree. If no
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
    warning("No sorted stable ID list found. Will return lexically sorted unique stable ID, which might not be what you wanted.");

    $self->_sort_stable_ids;
  }
  
  return $self->{'sorted_tree'}->{'stable_ids'};
}


=head2 _sort_stable_ids

  Example     : 
  Description : Returns a list of stable IDs in this history tree, sorted
                alphabetically. This is the simplest sort function used and
                doesn't try to untangle the tree.
  Return type : Arrayref
  Exceptions  : none
  Caller      : get_unique_stable_ids
  Status      : At Risk
              : under development

=cut

sub _sort_stable_ids {
  my $self = shift;

  unless ($self->{'sorted_tree'}->{'stable_ids'}) {
    my %unique = map { $_->stable_id => 1 } @{ $self->get_all_ArchiveStableIds };
    $self->{'sorted_tree'}->{'stable_ids'} = [sort keys %unique];
  }
  
  return $self->{'sorted_tree'}->{'stable_ids'};
}


=head2 coords_by_ArchiveStableId 

  Arg[1]      : Bio::EnsEMBL::ArchiveStableId $archive_id
                The ArchiveStableId to get tree grid coordinates for
  Example     : my ($x, $y) =
                  @{ $history_tree->coords_by_ArchiveStableId($archive_id) };
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
    unless ($archive_id->isa('Bio::EnsEMBL::ArchiveStableId'));
  
  return $self->{'sorted_tree'}->{'coords'}->{$self->_node_id($archive_id)}
    || [];
}


=head2 calculate_simple_coords 

  Example     : $history_tree->calculate_simple_coords;
  Description : Pre-calculates the grid coordinates of all nodes in the tree.
                This most simple method for this task sorts releases
                chronologically and stable IDs alphabetically. No efforts are
                made to disentangle the tree.
  Return type : none
  Exceptions  : none
  Caller      : ArchiveStableIdAdaptor::fetch_history_tree_by_stable_id
  Status      : At Risk
              : under development

=cut

sub calculate_simple_coords {
  my $self = shift;

  # reset any previous tree cordinate calculations
  $self->reset_tree;

  # the "master" information for the sorted tree is stored as the sorted lists
  # of releases (x) and stable IDs (y). Sort them now.
  my $db_names = $self->get_release_db_names;
  my $stable_ids = $self->_sort_stable_ids;
  
  # for performance reasons, additionally store coordinates in a lookup hash
  foreach my $archive_id (@{ $self->get_all_ArchiveStableIds }) {
  
    # coordinates are positions in the sorted lists
    my $x = $self->_index_of($archive_id->db_name, $db_names);
    my $y = $self->_index_of($archive_id->stable_id, $stable_ids);
  
    $self->{'sorted_tree'}->{'coords'}->{$self->_node_id($archive_id)} =
      [ $x, $y ];
  }
}


=head2 _index_of

  Arg[1]      : 
  Example     : my @array = (a, b, c);
                my $i = _index_of('b', \@array); # will return 1
  Description : Returns the index of an element in an array
  Return type : Int (or undef if element is not found in array)
  Exceptions  : thrown on wrong argument types
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub _index_of {
  my ($self, $element, $arrayref) = @_;

  throw("Expecting arrayref argument.") unless (ref($arrayref) eq 'ARRAY');

  my @array = @$arrayref;

  while (my $e = pop(@array)) {
    return scalar(@array) if ($e eq $element);
  }

  return undef;
}


=head2 reset_tree

  Example     : $history_tree->reset_tree;
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


1;


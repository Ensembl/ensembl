# EnsEMBL module for MarkerAdaptor
# Copyright EMBL-EBI/Sanger center 2002
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor

=head1 SYNOPSIS


=head1 DESCRIPTION

Provides database interaction for the Bio::EnsEMBL::Map::Marker object

=cut

package Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor;

use strict;

use vars ('@ISA');

use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MapLocation;
use Bio::EnsEMBL::Map::MarkerSynonym;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);




=head2 fetch_all

  Arg [1]    : none
  Example    : @all_markers = @{$marker_adaptor->fetch_all};
  Description: Retrieves all markers from the database
  Returntype : listref of Bio::EnsEMBL::Map::Markers
  Exceptions : none
  Caller     : general

=cut

sub fetch_all {
  my $self = shift;
  my $dbID = shift;

  my $sth = $self->prepare("SELECT m.marker_id, m.priority, m.left_primer, 
                                   m.right_primer, m.type,
                                   m.min_primer_dist, m.max_primer_dist,
                                   ms.marker_synonym_id, ms.name, ms.source
                            FROM   marker m 
                            LEFT JOIN marker_synonym ms 
                            ON     ms.marker_synonym_id = 
                                    m.display_marker_synonym_id");

  $sth->execute;

  my( $marker_id, $priority, $left_primer, $right_primer, $type,
      $min_pdist, $max_pdist, $ms_id, $ms_name, $ms_src);

  $sth->bind_columns(\$marker_id, \$priority,
		     \$left_primer, \$right_primer, \$type, \$min_pdist, \$max_pdist,
		     \$ms_id, \$ms_name, \$ms_src);

  my @out;
  while($sth->fetch) {
    #create a display marker synonym for each marker created, if one is defined
    my $synonym;
    if($ms_id) { 
      $synonym = Bio::EnsEMBL::Map::MarkerSynonym->new
	($ms_id, $ms_src, $ms_name);
    }	
    
    push @out, Bio::EnsEMBL::Map::Marker->new
      ($marker_id, $self, $left_primer, $right_primer, $min_pdist, $max_pdist,
       $priority, $type, $synonym);
  }

  return \@out;
}



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               The internal identifier of the Marker to retrieve
  Example    : $marker = $marker_adaptor->fetch_by_dbID(123);
  Description: Retrieves a marker object from the database via its internal
               identifier.
  Returntype : Bio::EnsEMBL::Map::Marker
  Exceptions : thrown if no marker with $dbID is present in the database
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  $self->throw('dbID argument is required') unless($dbID);

  my $sth = $self->prepare("SELECT m.marker_id, m.display_marker_synonym_id, 
                                   m.priority, m.left_primer, m.right_primer,
                                   m.type,
                                   m.min_primer_dist, m.max_primer_dist, 
                                   ms.marker_synonym_id,
                                   ms.source, ms.name
                            FROM marker m, marker_synonym ms
                            WHERE m.marker_id = ?
                            AND ms.marker_id = m.marker_id");

  $sth->execute($dbID);

  my( $marker_id, $display_ms_id, $priority, $left_primer, $right_primer,
      $type, $min_pdist, $max_pdist, $ms_id, $ms_src, $ms_name);

  $sth->bind_columns(\$marker_id, \$display_ms_id, \$priority, 
		     \$left_primer, \$right_primer, \$type, 
		     \$min_pdist, \$max_pdist, 
		     \$ms_id, \$ms_src, \$ms_name);

  my $display_synonym;
  my @synonyms;
  while($sth->fetch) {
    #create a new synonym for each row
    my $s = new Bio::EnsEMBL::Map::MarkerSynonym->new($ms_id, $ms_src, 
						     $ms_name);
    $display_synonym = $s if($display_ms_id == $ms_id);
    push @synonyms, $s;
  }

  $sth->finish;

  unless($marker_id) {
    $self->throw("marker with dbID=[$dbID] not present in database");
  }
    
  #now create the marker
  return new Bio::EnsEMBL::Map::Marker->new(
     $marker_id, $self, $left_primer, $right_primer,$min_pdist, $max_pdist,
     $priority, $type,$display_synonym, \@synonyms);
}


=head2 fetch_all_by_synonym

  Arg [1]    : string $synonym
               An name of this marker
  Arg [2]    : (opional) string $source
               The source of this name
  Example    : @markers = @{$marker_adaptor->fetch_all_by_synonym($id)};
  Description: Retrieves a list of markers with the synonym (alias) $synonym
               and from source $source.  In most cases the list will have a 
               single element, however it is possible that multiple markers 
               with the same synonym exist.
  Returntype : listref of Bio::EnsEMBL::Map::Markers
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_synonym {
  my ($self, $synonym, $source) = @_;

  $self->throw("synonym argument is required") unless($synonym);

  my $q = "SELECT marker_id 
           FROM   marker_synonym ms
           WHERE  ms.name = ?";

  my @bind_vals = ($synonym);

  if($source) {
    $q .= " AND ms.source = ?";
    push(@bind_vals, $source);
  }

  my $sth = $self->prepare($q);
  $sth->execute(@bind_vals);

  my @out = ();
  my $marker_id;
  my %seen;
  
  #fetch the markers and filter out duplictes
  while(($marker_id) = $sth->fetchrow_array) {
    next if $seen{$marker_id};

    push @out, $self->fetch_by_dbID($marker_id);
    $seen{$marker_id} = 1;
  }
  
  $sth->finish;

  return \@out;
}
  


=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::Map::Marker $marker
  Example    : $marker_adaptor->fetch_attributes($marker);
  Description: Fetches the marker_synonym and map_location attributes of
               a marker.  This is done so that these attributes can be
               lazy-loaded on request. 
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::Map::Marker::marker

=cut

sub fetch_attributes {
  my $self = shift;
  my $marker = shift;

  my $m_id = $marker->dbID;

  $self->throw('Marker argument does not have a dbID') unless($m_id);

  #
  # First Retrieve synonyms for this marker
  #
  $marker->flush_MarkerSynonyms;

  my $sth = $self->prepare("SELECT ms.marker_synonym_id, ms.source, ms.name
                            FROM   marker_synonym ms
                            WHERE  ms.marker_id = ?");

  my @syns = ();
  my ($ms_id, $ms_src, $ms_name);

  $sth->execute($m_id);
  $sth->bind_columns(\$ms_id, \$ms_src, \$ms_name);

  while($sth->fetch) {
    push @syns, Bio::EnsEMBL::Map::MarkerSynonym->new($ms_id,$ms_src,$ms_name);
  }
  $sth->finish;

  $marker->add_MarkerSynonyms(@syns) if(@syns);

  #
  # Now retrieve map locations for this marker
  #
  $marker->flush_MapLocations;

  $sth = $self->prepare("SELECT mloc.chromosome_name, mloc.position,
                                mloc.lod_score, map.map_name, ms.name
                         FROM   marker_map_location mloc, map map,
                                marker_synonym ms
                         WHERE  mloc.marker_id = ? 
                         AND    map.map_id = mloc.map_id
                         AND   ms.marker_synonym_id = mloc.marker_synonym_id");

  my($chr_name, $pos, $lod, $mname, $name);
  my @mlocs;

  $sth->execute($m_id);

  $sth->bind_columns(\$chr_name, \$pos, \$lod, \$mname, \$name);

  while($sth->fetch) {
    push(@mlocs, Bio::EnsEMBL::Map::MapLocation->new($name, $mname,
                                                     $chr_name,$pos,$lod));
  }

  $sth->finish;

  $marker->add_MapLocations(@mlocs);
}


1;

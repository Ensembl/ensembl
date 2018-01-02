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

Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor

=head1 SYNOPSIS

=head1 DESCRIPTION

Provides database interaction for the Bio::EnsEMBL::Map::Marker object

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor;

use strict;

use vars ('@ISA');

use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MapLocation;
use Bio::EnsEMBL::Map::MarkerSynonym;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);




=head2 fetch_all

  Arg [1]    : none
  Example    : @all_markers = @{$marker_adaptor->fetch_all};
  Description: Retrieves all markers from the database
  Returntype : listref of Bio::EnsEMBL::Map::Markers
  Exceptions : none
  Caller     : general
  Status     : stable

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
  Status     : stable

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
    $self-> warning("marker with dbID=[$dbID] not present in database");
    return undef;
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
  Status     : stable

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

    # some synonyms point to markers that don't exist, so only add genuine ones
    my $marker = $self->fetch_by_dbID($marker_id);
    push @out, $marker if ($marker);
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
  Status     : stable

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



=head2 store

  Arg [1]    : Bio::EnsEMBL::Map::Marker
  Example    : $marker_adaptor->store(@markers);
  Description: Stores a list of markers in this database.
               The dbID and adaptor of each marker will be set on successful
               storing.
  Returntype : 1 on success
  Exceptions : thrown if not all data needed for storing is populated in the
               marker
  Caller     : general
  Status     : stable

=cut

sub store {
  my ($self, @markers) = @_;

  MARKER:foreach my $marker( @markers ){

    if($marker->dbID){
      if($self->fetch_by_dbID($marker->dbID)){
        next MARKER;
      }
    }
    #
    # Sanity check
    #
    if(!$marker ||
       !ref($marker) || 
       !$marker->isa('Bio::EnsEMBL::Map::Marker')) {
      $self->throw('Incorrect argument [$marker] to store.  Expected ' .
                   'Bio::EnsEMBL::Map::Marker');
    }

    # Don't store if already stored
    if($marker->is_stored($self->db())) {
      warning('Marker ['.$marker->dbID.'] is already stored in this DB.');
      next;
    }

    # Get/test the display marker synonym
    my $display_synonym = $marker->display_MarkerSynonym;
    if(!$display_synonym || !ref($display_synonym) ||
       !$display_synonym->isa('Bio::EnsEMBL::Map::MarkerSynonym')) {
      $self->throw('Cannot store Marker without an associated '.
                   'display_MarkerSynonym');
    }

    # Store the Marker itself
    my $q = qq(
INSERT INTO marker ( left_primer, right_primer,
                     min_primer_dist, max_primer_dist, 
                     priority, type)
            VALUES ( ?,?,?,?,?,?) );  

    my $sth = $self->prepare($q);

    $sth->execute( $marker->left_primer      || '',
                   $marker->right_primer     || '',
                   $marker->min_primer_dist  || 0,
                   $marker->max_primer_dist  || 0,
                   $marker->priority,
                   $marker->type );

    my $dbID = $self->last_insert_id('marker_id', undef, 'marker');
    $marker->dbID($dbID);
    $marker->adaptor($self);

    if(!$display_synonym->dbID) {
      # Store synonym
      $self->_store_MarkerSynonym($marker,$display_synonym);
    }
    my $display_synonym_id = $display_synonym->dbID ||
        $self->throw('display_MarkerSynonym must have dbID to store Marker');

    # Update the marker with the display synonym
    my $qup = qq(
UPDATE marker 
SET    display_marker_synonym_id = $display_synonym_id
                 WHERE  marker_id = ? );
    my $sthup = $self->prepare($qup);
    $sthup->execute($dbID);

    # Loop through all MarkerSynonyms and store if needed
    foreach my $synonym( @{$marker->get_all_MarkerSynonyms} ){
      if(!$synonym->dbID) {
        $self->_store_MarkerSynonym($marker,$synonym);
      }
    }

    # Loop through all MapLocations and store if needed
    foreach my $loc( @{$marker->get_all_MapLocations} ){
      # Dunno how to implement this :( Just bomb out
      $self->throw( 'Storing of MapLocation objects is not yet implemented' )
    }

  }
  return 1;
}


=head2 _store_MarkerSynonym

  Arg [1]    : Bio::EnsEMBL::Map::Marker
  Arg [2]    : Bio::EnsEMBL::Map::MarkerSynonym
  Example    : $marker_adaptor->_store_MarkerSynonym($marker,$ms);
  Description: Stores a marker synonym attached to the marker in the database
               The dbID of each MarkerSynonym will be set on successful
               storing.
  Returntype : dbID of the MarkerSynonym
  Exceptions : thrown if not all data needed for storing is populated
  Caller     : $self->store
  Status     : stable

=cut

sub _store_MarkerSynonym{
  my $self   = shift;
  my $marker = shift;
  my $ms     = shift;

  # Sanity check
  if(!$marker || !ref($marker) || 
     !$marker->isa('Bio::EnsEMBL::Map::Marker')) {
    $self->throw("Incorrect argument [$marker] to _store_MarkerSynonym." .
                 "Expected Bio::EnsEMBL::Map::Marker");
  }
  if(!$ms || !ref($ms) || 
     !$ms->isa('Bio::EnsEMBL::Map::MarkerSynonym')) {
    $self->throw("Incorrect argument [$ms] to _store_MarkerSynonym." .
                 "Expected Bio::EnsEMBL::Map::MarkerSynonym");
  }

  # Don't store if already stored
  if($ms->dbID) {
    warning('MarkerSynonym ['.$ms->dbID.'] is already stored in this DB.');
    return;
  }

  my $marker_id = $marker->dbID ||
      throw( "Marker has no dbID. Cannot store MarkerSynonym" );

  # Store the synonym
  my $q = qq(
INSERT INTO marker_synonym ( marker_id, source, name )
VALUES ( ?,?,?) );

  my $sth = $self->prepare($q);

  $sth->execute( $marker_id,
                 $ms->source,
                 $ms->name );

  my $dbID = $self->last_insert_id('marker_synonym_id', undef, 'marker_synonym');
  $ms->dbID($dbID);
  return $dbID;
}

1;

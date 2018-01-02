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

Bio::EnsEMBL::Map::Ditag

=head1 SYNOPSIS

  my $feature = Bio::EnsEMBL::Map::Ditag->new(
    -dbID      => $tag_id,
    -name      => $name,
    -type      => $type,
    -tag_count => $tag_count,
    -sequence  => $sequence,
    -adaptor   => $dbAdaptor
  );

=head1 DESCRIPTION

Represents an unmapped ditag object in the EnsEMBL database.
Corresponds to original tag containing the full sequence. This can be
a single piece of sequence like CAGE tags or a ditag with concatenated
sequence from 5' and 3' end like GIS or GSC tags.

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::Ditag;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [1]    : (optional) int $dbID
  Arg [2]    : (optional) string name
  Arg [3]    : (optional) string type
  Arg [4]    : (optional) int tag_count
  Arg [5]    : (optional) string sequence
  Arg [6]    : (optional) Bio::EnsEMBL::Map::DBSQL::DitagAdaptor $adaptor

  Description: Creates a new ditag
  Returntype : Bio::EnsEMBL::Map::Ditag
  Exceptions : none
  Caller     : general

=cut

sub new {
  my ($caller, @args) = @_;
  my ($dbID, $name, $type, $tag_count, $sequence, $adaptor) = rearrange(
      [ 'DBID', 'NAME', 'TYPE', 'TAG_COUNT', 'SEQUENCE', 'ADAPTOR' ], @args);
  my $class = ref($caller) || $caller;

  if(!$name or !$type or !$sequence) {
    throw('Missing information for Ditag object:
              Bio::EnsEMBL::Map::Ditag->new (
                                              -dbID      => $tag_id,
                                              -name      => $name,
                                              -type      => $type,
                                              -tag_count => $tag_count,
                                              -sequence  => $sequence,
                                              -adaptor   => $dbAdaptor
                                             );');
  }

  if(!$tag_count){ $tag_count = 0; }

  if(!($sequence =~ /^[ATCGN]+$/i)){
    throw('ditag sequence contains non-standard characters: '.$sequence);
  }

  my $self = bless( {'dbID'        => $dbID,
                     'name'        => $name,
                     'type'        => $type,
		     'tag_count'   => $tag_count,
		     'sequence'    => $sequence
                    }, $class);

  $self->adaptor($adaptor);
  return $self;
}

=head2 name

  Arg [1]    : (optional) string $type
  Example    : $type = $ditag->name;
  Description: Getter/Setter for the name of a ditag
  Returntype : text
  Caller     : general

=cut

sub name {
  my $self = shift;

  if(@_) {
    $self->{'name'} = shift;
  }

  return $self->{'name'};
}

=head2 dbID

  Arg [1]    : (optional) int id
  Example    : $ditag_id = $ditag->dbID;
  Description: Getter/Setter for the dbID of a ditag
  Returntype : int
  Caller     : general

=cut

sub dbID {
  my $self = shift;

  if(@_) {
    $self->{'dbID'} = shift;
  }

  return $self->{'dbID'};
}


=head2 type

  Arg [1]    : (optional) string $type
  Example    : $type = $ditag->type;
  Description: Getter/Setter for the type of a ditag
  Returntype : text
  Caller     : general

=cut

sub type {
  my $self = shift;

  if(@_) {
    $self->{'type'} = shift;
  }

  return $self->{'type'};
}

=head2 tag_count

  Arg [1]    : (optional) string $tag_count
  Example    : $type = $ditag->tag_count;
  Description: Getter/Setter for the tag_count of a ditag
  Returntype : int
  Caller     : general

=cut

sub tag_count {
  my $self = shift;

  if(@_) {
    $self->{'tag_count'} = shift;
  }

  return $self->{'tag_count'};
}

=head2 sequence

  Arg [1]    : (optional) string $sequence
  Example    : $sequence = $ditag->sequence;
  Description: Getter/Setter for the sequence of a ditag
  Returntype : text
  Caller     : general

=cut

sub sequence {
  my $self = shift;

  if(@_) {
    $self->{'sequence'} = shift;
  }

  return $self->{'sequence'};
}


=head2 get_ditagFeatures

  Arg        : none
  Example    : @features = @{$ditag->get_ditagFeatures};
  Description: Fetch ditag_features created from this ditag
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeature
  Caller     : general

=cut

sub get_ditagFeatures {
  my $self = shift;

  return $self->adaptor->db->get_adaptor("ditagFeature")
          ->fetch_all_by_ditagID($self->dbID);
}

1;

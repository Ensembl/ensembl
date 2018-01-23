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

Bio::EnsEMBL::Map::DBSQL::DitagAdaptor

=head1 SYNOPSIS

  my $ditagadaptor = $db->get_DitagAdaptor();
  my @ditags       = @{ $ditagadaptor->fetch_by_type("ZZ11") };

=head1 DESCRIPTION

Provides database interaction for the Bio::EnsEMBL::Map::Ditag object

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::DBSQL::DitagAdaptor;

use strict;
use vars ('@ISA');

use Bio::EnsEMBL::Map::Ditag;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_all_by_name

  Arg [1]    : ditag name
  Example    : $tag = $ditag_adaptor->fetch_by_name("U3");
  Description: Retrieves ditags from the database using the name
  Returntype : listref of Bio::EnsEMBL::Map::Ditag
  Caller     : general

=cut

sub fetch_all_by_name {
  my ($self, $tagname) = @_;

  if(!$tagname){
    throw "must be called with a name of a ditag.";
  }
  my $sth = $self->prepare("SELECT d.ditag_id, d.name, d.type, d.tag_count, d.sequence
                            FROM   ditag d
                            WHERE  d.name = ?");
  $sth->execute($tagname);
  my $result = $self->_fetch($sth);

  return $result;

}


=head2 fetch_by_dbID

  Arg [1]    : ditag type
  Example    : @all_tags = @{$ditag_adaptor->fetch_by_dbID(1003)};
  Description: Retrieve ditags with a certain id from the database
  Returntype : Bio::EnsEMBL::Map::Ditag
  Caller     : general

=cut

sub fetch_by_dbID {
  my ($self, $tagid) = @_;

  if(!$tagid){
    throw "must be called with the type of a ditag.";
  }
  my $sth = $self->prepare("SELECT d.ditag_id, d.name, d.type, d.tag_count, d.sequence
                            FROM   ditag d
                            WHERE  d.ditag_id = ?");
  $sth->execute($tagid);
  my $result = $self->_fetch($sth);

  return $result->[0];
}


=head2 fetch_all_by_type

  Arg [1]    : ditag type
  Example    : @all_tags = @{$ditag_adaptor->fetch_by_type("SME005")};
  Description: Retrieves all ditags of a certain type from the database
  Returntype : listref of Bio::EnsEMBL::Map::Ditag
  Caller     : general

=cut

sub fetch_all_by_type {
  my ($self, $tagtype) = @_;

  if(!$tagtype){
    throw "must be called with the type of a ditag.";
  }
  my $sth = $self->prepare("SELECT d.ditag_id, d.name, d.type, d.tag_count, d.sequence
                            FROM   ditag d
                            WHERE  d.type = ?");
  $sth->execute($tagtype);
  my $result = $self->_fetch($sth);

  return $result;
}


=head2 fetch_by_name_and_type

  Arg [1]    : ditag name
  Arg [2]    : ditag type
  Example    : $tag = $ditag_adaptor->fetch_by_name_and_type("U3", "SME005");
  Description: Retrieves ditags from the database using name/type combination
               which should be non-ambiguous
  Returntype : Bio::EnsEMBL::Map::Ditag
  Caller     : general

=cut

sub fetch_by_name_and_type {
  my ($self, $tagname, $tagtype) = @_;

  if(!$tagname or !$tagtype){
    throw "must be called with a name and type of a ditag.";
  }
  my $sth = $self->prepare("SELECT d.ditag_id, d.name, d.type, d.tag_count, d.sequence
                            FROM   ditag d
                            WHERE  d.name = ? and d.type = ?");
  $sth->execute($tagname, $tagtype);
  my $result = $self->_fetch($sth);

  return $result->[0];
}


=head2 fetch_all

  Args       : none
  Example    : @all_tags = @{$ditag_adaptor->fetch_all};
  Description: Retrieves all ditags from the database
  Returntype : listref of Bio::EnsEMBL::Map::Ditag
  Caller     : general

=cut

sub fetch_all {
  my ($self) = @_;

  my $sth = $self->prepare("SELECT d.ditag_id, d.name, d.type, d.tag_count, d.sequence
                            FROM   ditag d");
  $sth->execute;
  my $result = $self->_fetch($sth);

  return $result;
}


=head2 fetch_all_with_limit

  Arg [1]    : ditag type
  Arg [2]    : row limit
  Arg [3]    : row offset
  Description: fetch_by_type with row limit and offset
  Returntype : listref of Bio::EnsEMBL::Map::Ditag
  Caller     : general

=cut

sub fetch_all_with_limit {
  my ($self, $tagtype, $limit, $offset) = @_;

  my @ditags = ();
  my $sql = "SELECT d.ditag_id, d.name, d.type, d.tag_count, d.sequence ".
            "FROM ditag d ".
            "WHERE d.type = ? LIMIT ? OFFSET ?;";
  my $sth = $self->db->dbc->prepare($sql);
  $sth->execute($tagtype, $limit, $offset);
  my $result = $self->_fetch($sth);

  return $result;
}


=head2 _fetch

  Arg [1]    : statement handler object
  Description: generic sql-fetch function for the ditag fetch methods
  Returntype : listref of Bio::EnsEMBL::Map::Ditag
  Caller     : private

=cut

sub _fetch {
  my ($self, $sth) = @_;

  my($tag_id, $name, $type, $count, $sequence);
  my @tags;

  $sth->bind_columns(\$tag_id, \$name, \$type, \$count, \$sequence);
  while($sth->fetch) {
    push @tags, Bio::EnsEMBL::Map::Ditag->new (
                          -dbID      => $tag_id,
                          -name      => $name,
                          -type      => $type,
		          -tag_count => $count,
                          -sequence  => $sequence,
                          -adaptor   => $self,
                          );
  }

  return \@tags;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::Map::Ditag
  Example    : $ditag_adaptor->store(\@ditags);
  Description: Stores a single ditag or
               a list of ditags in this database.
  Returntype : none
  Caller     : general

=cut

sub store {
  my ($self, $ditags) = @_;

  if(ref $ditags eq 'ARRAY'){
    if(scalar(@$ditags) == 0){
      throw("Must call store with ditag or list ref of ditags");
    }
  }
  elsif($ditags){
    my @ditags;
    push @ditags, $ditags;
    $ditags = \@ditags;
  }
  else{
    throw("Must call store with ditag or list ref of ditags not ".$ditags);
  }

  my $db = $self->db() or throw "Couldn t get database connection.";

 TAG:
  foreach my $ditag (@$ditags) {

    if ( !ref $ditag || !$ditag->isa("Bio::EnsEMBL::Map::Ditag") ) {
      throw( "Object must be an Ensembl Ditag, " . "not a [" . ref($ditag) . "]" );
    }

    if ( $ditag->is_stored($db) ) {
      warning( "Ditag [" . $ditag->dbID . "] is already stored in this database." );
      next TAG;
    }

    #check if tag with same name/type already exists
    my $sth = $self->prepare( "SELECT COUNT(*) FROM ditag 
                               WHERE name = ? AND type = ?" );
    $sth->execute($ditag->name, $ditag->type);
    if($sth->fetchrow() > 0){
      warning( "Ditag with name/type ".$ditag->name." / ".$ditag->type.
               " is already stored in this database.\n".
	       "Use update_ditag() instead.");
      next TAG;
    }

    if ( $ditag->dbID ) {
      my $sth = $self->prepare( "INSERT INTO ditag( ditag_id , name, type, tag_count, sequence ) ".
			        "VALUES( ?,?,?,?,? )" );
      $sth->bind_param(1,$ditag->dbID,SQL_INTEGER);
      $sth->bind_param(2,$ditag->name,SQL_VARCHAR);
      $sth->bind_param(3,$ditag->type,SQL_VARCHAR);
      $sth->bind_param(4,$ditag->tag_count,SQL_VARCHAR);
      $sth->bind_param(5,$ditag->sequence,SQL_VARCHAR);
      $sth->execute();
    } else {
      my $sth = $self->prepare( "INSERT INTO ditag( name, type, tag_count, sequence ) ".
			        "VALUES( ?,?,?,? )" );
      $sth->bind_param(1,$ditag->name,SQL_VARCHAR);
      $sth->bind_param(2,$ditag->type,SQL_VARCHAR);
      $sth->bind_param(3,$ditag->tag_count,SQL_VARCHAR);
      $sth->bind_param(4,$ditag->sequence,SQL_VARCHAR);
      $sth->execute();
      my $dbID = $self->last_insert_id('ditag_id', undef, 'ditag');
      $ditag->dbID($dbID);
      $ditag->adaptor($self);
    }
  }

  return 1;
}


=head2 print_creation

  Arg [1]    : ditag probe name
  Arg [2]    : ditag type
  Arg [3]    : ditag count
  Arg [4]    : ditag sequence
  Arg [5]    : (optional) ditag dbID
  Description: convenience method to produce SQL insert statements
               to speed up the creation of new ditags
  Returntype : string

=cut

sub print_creation {
  my ($self, $probe_name, $type, $count, $sequence, $dbid) = @_;
  my $string;
  if($dbid){
    $string = "INSERT INTO ditag( ditag_id, name, type, tag_count, sequence ) ".
	      "VALUES($dbid, '".$probe_name."', '".$type."', ".$count."'".$sequence."');\n";
  }
  else {
    $string = "INSERT INTO ditag( name, type, tag_count, sequence ) ".
	      "VALUES('".$probe_name."', '".$type."', ".$count.", '".$sequence."');\n";
  }

  return $string;
}


=head2 update_ditag

  Arg [1]    : ditag to update
  Description: update an existing ditag with new values
  Returntype : true on success

=cut

sub update_ditag {
  my ($self, $ditag) = @_;

  my $sth = $self->prepare( "UPDATE ditag SET name=?, type=?, tag_count=?, sequence=? where ditag_id=?;" );
  my $result =$sth->execute(
                            $ditag->name,
                            $ditag->type,
                            $ditag->tag_count,
                            $ditag->sequence,
                            $ditag->dbID,
                           );

  return $result;
}


=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$da->list_dbIDs()};
  Description: Gets an array of internal IDs for all Ditag objects in
               the current database.
  Returntype : List of ints
  Exceptions : None

=cut

sub list_dbIDs {
  my ($self, $ordered) = @_;
	
  return $self->_list_dbIDs('ditag');
}

1;

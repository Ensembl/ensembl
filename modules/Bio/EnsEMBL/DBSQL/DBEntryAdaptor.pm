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

Bio::EnsEMBL::DBSQL::DBEntryAdaptor -
MySQL Database queries to load and store external object references.

=head1 SYNOPSIS

  $db_entry_adaptor =
    $registry->get_adaptor( 'Human', 'Core', 'DBEntry' );

  $db_entry = $db_entry_adaptor->fetch_by_dbID($id);

  my $gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );

  my $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000101367');

  @db_entries = @{ $db_entry_adaptor->fetch_all_by_Gene($gene) };
  @gene_ids   = $db_entry_adaptor->list_gene_ids_by_extids('BAB15482');

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::DBEntryAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::IdentityXref;
use Bio::EnsEMBL::OntologyXref;

use Bio::EnsEMBL::Utils::Exception qw(deprecate throw warning);

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               the unique database identifier for the DBEntry to retrieve
  Example    : my $db_entry = $db_entry_adaptor->fetch_by_dbID($dbID);
  Description: Retrieves a dbEntry from the database via its unique
               identifier.
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my ( $self, $dbID ) = @_;

  my $sth = $self->prepare(
    "SELECT  xref.xref_id,
            xref.dbprimary_acc,
            xref.display_label,
            xref.version,
            exDB.priority,
            exDB.db_name,
            exDB.db_display_name,
            exDB.db_release,
            es.synonym,
            xref.info_type,
            xref.info_text,
            exDB.type,
            exDB.secondary_db_name,
            exDB.secondary_db_table,
            xref.description
    FROM    (xref, external_db exDB)
    LEFT JOIN external_synonym es ON
            es.xref_id = xref.xref_id
    WHERE   xref.xref_id = ?
    AND     xref.external_db_id = exDB.external_db_id" );

  $sth->bind_param( 1, $dbID, SQL_INTEGER );
  $sth->execute();

  my $exDB;

  my $max_rows = 1000;

  while ( my $rowcache = $sth->fetchall_arrayref( undef, $max_rows ) ) {
      #$description refers to the external_db description, while $desc was referring the xref description
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $refID,               $dbprimaryId,
           $displayid,           $version,
           $priority,
           $dbname,              $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      if ( !defined($exDB) ) {
        $exDB =
          Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $dbname,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
			               -description         => $description
          );


      }

      if ( defined($synonym) ) { $exDB->add_synonym($synonym) }

    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return $exDB;
} ## end sub fetch_by_dbID


sub _get_all_dm_loc_sth {
  my ($self, $constraint ,$ensembl_object ) = @_;
  my $object_type;
  if($ensembl_object->isa("Bio::EnsEMBL::Gene")){
    $object_type = "Gene";
  }
  elsif($ensembl_object->isa("Bio::EnsEMBL::Transcript")){
    $object_type = "Transcript";
  }
  elsif($ensembl_object->isa("Bio::EnsEMBL::Translation")){
    $object_type = "Translation";
  }
  elsif($ensembl_object->isa("Bio::EnsEMBL::Operon")){
    $object_type = "Operon";
  }
  elsif($ensembl_object->isa("Bio::EnsEMBL::OperonTranscript")){
    $object_type = "OperonTranscript";
  }
  else{
    warn(ref($ensembl_object)." is not a Gene Transcript or Translation object??\n");
    return undef;
  }
  my $sql = "SELECT xref.xref_id,
               xref.dbprimary_acc,
            xref.display_label,
            xref.version,
            exDB.priority,
            exDB.db_name,
            exDB.db_display_name,
            exDB.db_release,
            es.synonym,
            xref.info_type,
            xref.info_text,
            exDB.type,
            exDB.secondary_db_name,
            exDB.secondary_db_table,
            xref.description
    FROM    (xref, external_db exDB, dependent_xref dx, object_xref ox)
    LEFT JOIN external_synonym es ON
            es.xref_id = xref.xref_id
    WHERE   xref.external_db_id = exDB.external_db_id AND
            ox.xref_id = xref.xref_id AND
            ox.ensembl_object_type = \'$object_type\' AND
            ox.ensembl_id = ".$ensembl_object->dbID();

  if($constraint){
    $sql .= " AND $constraint";
  }
  else{
    die "NO constraint???\n";
  }

  my $sth = $self->prepare($sql) || die "Could not prepare $sql";

  return $self->_get_all_dm($sth);
}

sub _get_all_dm_sth {
  my ( $self, $constraint) = @_;

 my $sql = "SELECT xref.xref_id,
               xref.dbprimary_acc,
            xref.display_label,
            xref.version,
            exDB.priority,
            exDB.db_name,
            exDB.db_display_name,
            exDB.db_release,
            es.synonym,
            xref.info_type,
            xref.info_text,
            exDB.type,
            exDB.secondary_db_name,
            exDB.secondary_db_table,
            xref.description
    FROM    (xref, external_db exDB, dependent_xref dx)
    LEFT JOIN external_synonym es ON
            es.xref_id = xref.xref_id
    WHERE   xref.external_db_id = exDB.external_db_id ";

  if($constraint){
    $sql .= "AND $constraint";
  }
  else{
    die "NO constraint???\n";
  }

  my $sth = $self->prepare($sql) || die "Could not prepare $sql";

  return $self->_get_all_dm($sth);
}


sub _get_all_dm{

  my ($self, $sth) = @_;

#  $sth->bind_param( 1, $dm_dbid, SQL_INTEGER );

#  print $sth."\n";
  $sth->execute() || die "Not able to execute statement handle";

  my @list =();
  my %seen;

  my $max_rows = 1000;
  while ( my $rowcache = $sth->fetchall_arrayref(undef, $max_rows) ) {
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $dbID,                $dbprimaryId,
           $displayid,           $version,
           $priority,
           $dbname,              $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      if ( !defined($seen{$dbID}) ) {
       my $exDB =
          Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $dbname,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
			   -description         => $description
          );

	if ($synonym) { $exDB->add_synonym($synonym) };
	$seen{$dbID} = 1;
	push @list, $exDB;
      }



    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return \@list;

}


=head2 get_all_dependents

  Args[1]    : dbID of the DBentry to get the dependents of.
  Args[2]    : (optional) Bio::EnsEMBL::Gene, Transcript or Translation object
  Example    : my @dependents = @{ $dbe_adaptor->get_all_dependents(1234) };
  Description: Get a list of DBEntrys that are depenednet on the DBEntry.
               if an ensembl gene transcript or translation is given then only
               the ones on that object will be given
  Returntype : listref of DBEntrys. May be empty.
  Exceptions : none
  Caller     : DBEntry->get_all_dependnets
  Status     : UnStable

=cut

sub get_all_dependents {
  my ( $self, $dbid, $ensembl_object) = @_;
  
  if(defined($ensembl_object) and !($ensembl_object->isa("Bio::EnsEMBL::Feature") or $ensembl_object->isa("Bio::EnsEMBL::Translation"))){
    die ref($ensembl_object)." is not an Gene Transcript or Translation";
  }
  
  my $constraint = " dx.master_xref_id = $dbid AND  dx.dependent_xref_id = xref.xref_id";
  if(defined($ensembl_object)){
    return $self->_get_all_dm_loc_sth($constraint, $ensembl_object);
  }
  else{
    return $self->_get_all_dm_sth($constraint, $ensembl_object);
  }

}

=head2 get_all_masters

  Args[1]    : dbID of the DBentry to get the masters of.
  Args[2]    : (optional) Bio::EnsEMBL::Gene, Transcript or Translation object
  Example    : my @masters = @{ $dbe_adaptor->get_all_masters(1234) };
  Description: Get a list of DBEntrys that are the masters of the DBEntry.
               if an ensembl gene transcript or translation is given then only
               the ones on that object will be given.
  Returntype : listref of DBEntrys. May be empty.
  Exceptions : none
  Caller     : DBEntry->get_all_masters
  Status     : UnStable

=cut

sub get_all_masters {
  my ( $self, $dbid, $ensembl_object ) = @_;
  
  if(defined($ensembl_object) and !($ensembl_object->isa("Bio::EnsEMBL::Feature") or $ensembl_object->isa("Bio::EnsEMBL::Translation"))){
    die ref($ensembl_object)." is not an Gene Transcript or Translation";
  }
  
  my $constraint = "dx.dependent_xref_id = $dbid AND  dx.master_xref_id = xref.xref_id";

  if(defined($ensembl_object)){
    return $self->_get_all_dm_loc_sth($constraint, $ensembl_object);
  }
  else{
    return $self->_get_all_dm_sth($constraint, $ensembl_object);
  }
#  return $self->_get_all_dm($constraint, $ensembl_object);
}


=head fetch_all_by_name

  Arg [1]    : string $name - The name of the external reference.
               found in accession, display_label or synonym
  Arg [2]    : (optional) string $dbname  - The name of the database which 
               the provided name is for.

  Example    : my $xref = @{$dbea->fetch_all_by_name('BRAC2','HGNC')}[0];
               print $xref->description(), "\n" if($xref);
  Description: Retrieves list of DBEntrys (xrefs) via a name.
               The accesion is looked for first then the synonym and finally
               the display_label.
               NOTE $dbname this is optional but adding this speeds the
               process up if you know what you are looking for.

               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : Bio::EnsEMBL::DBSQL::DBEntry
  Exceptions : thrown if arguments are incorrect
  Caller     : general, domainview
  Status     : Stable

=cut

sub fetch_all_by_name {
  my ( $self, $name, $dbname ) = @_;

  my $sql = (<<SQL);
  SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label, xref.version,
         exDB.priority, exDB.db_name, exDB.db_display_name, exDB.db_release,
            es.synonym, xref.info_type, xref.info_text,
            exDB.type, exDB.secondary_db_name, exDB.secondary_db_table,
            xref.description
    FROM    (xref, external_db exDB)
    LEFT JOIN external_synonym es ON
            es.xref_id = xref.xref_id
    WHERE  (xref.dbprimary_acc = ? or xref.display_label = ?)
    AND    xref.external_db_id = exDB.external_db_id
SQL

  if(defined $dbname){
    $sql .= " AND    exDB.db_name = ?";
  }
  my $sth = $self->prepare($sql);
  $sth->bind_param( 1, $name, SQL_VARCHAR );
  $sth->bind_param( 2, $name, SQL_VARCHAR );
  if(defined $dbname){
    $sth->bind_param( 3 , $dbname,    SQL_VARCHAR );
  }
  $sth->execute();

  my $max_rows = 1000;

  my $precache = $sth->fetchall_arrayref( undef, $max_rows ); # need to fetch to ensure rows() works for SQLite

  if ( !$sth->rows() && lc($dbname) eq 'interpro' ) {
  # This is a minor hack that means that results still come back even
  # when a mistake was made and no interpro accessions were loaded into
  # the xref table.  This has happened in the past and had the result of
  # breaking domainview

    $precache = undef;
    $sth->finish();
    $sth = $self->prepare(
      "SELECT   NULL,
                i.interpro_ac,
                i.id,
                NULL,
                NULL,
                'Interpro',
                NULL,
                NULL
        FROM    interpro i
        WHERE   i.interpro_ac = ?" );

    $sth->bind_param( 1, $name, SQL_VARCHAR );
    $sth->execute();
  }

  my %exDB;
  my @exDBlist;

  while ( my $rowcache = $precache || $sth->fetchall_arrayref( undef, $max_rows ) ) {
    $precache = undef;
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $dbID,                $dbprimaryId,
           $displayid,           $version,
           $priority,
           $dbname,              $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      if ( !defined $exDB{$dbID} ) {
	my $entrie = 
          Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $dbname,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
			   -description         => $description
          );
	$exDB{$dbID} = $entrie;
	push @exDBlist, $entrie;
      }
      if ($synonym) { $exDB{$dbID}->add_synonym($synonym) }

    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return \@exDBlist;
} ## end sub fetch_all_by_name



=head2 fetch_by_db_accession

  Arg [1]    : string $dbname - The name of the database which the provided
               accession is for.
  Arg [2]    : string $accession - The accesion of the external reference to
               retrieve.
  Example    : my $xref = $dbea->fetch_by_db_accession('Interpro','IPR003439');
               print $xref->description(), "\n" if($xref);
  Description: Retrieves a DBEntry (xref) via the name of the database
               it is from and its primary accession in that database.
               Undef is returned if the xref cannot be found in the
               database.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : Bio::EnsEMBL::DBSQL::DBEntry
  Exceptions : thrown if arguments are incorrect
  Caller     : general, domainview
  Status     : Stable

=cut

sub fetch_by_db_accession {
  my ( $self, $dbname, $accession ) = @_;

  my $sth = $self->prepare(
    "SELECT xref.xref_id,
            xref.dbprimary_acc,
            xref.display_label,
            xref.version,
            exDB.priority,
            exDB.db_name,
            exDB.db_display_name,
            exDB.db_release,
            es.synonym,
            xref.info_type,
            xref.info_text,
            exDB.type,
            exDB.secondary_db_name,
            exDB.secondary_db_table,
            xref.description
    FROM    (xref, external_db exDB)
    LEFT JOIN external_synonym es ON
            es.xref_id = xref.xref_id
    WHERE  xref.dbprimary_acc = ?
    AND    exDB.db_name like ?
    AND    xref.external_db_id = exDB.external_db_id" );

  $sth->bind_param( 1, $accession, SQL_VARCHAR );
  $sth->bind_param( 2, $dbname,    SQL_VARCHAR );
  $sth->execute();

  my $max_rows = 1000;

  my $precache = $sth->fetchall_arrayref( undef, $max_rows ); # need to fetch to ensure rows() works for SQLite

  if ( !$sth->rows() && lc($dbname) eq 'interpro' ) {
  # This is a minor hack that means that results still come back even
  # when a mistake was made and no interpro accessions were loaded into
  # the xref table.  This has happened in the past and had the result of
  # breaking domainview

    $precache = undef;
    $sth->finish();
    $sth = $self->prepare(
      "SELECT   NULL,
                i.interpro_ac,
                i.id,
                NULL,
                NULL,
                'Interpro',
                NULL,
                NULL
        FROM    interpro i
        WHERE   i.interpro_ac = ?" );

    $sth->bind_param( 1, $accession, SQL_VARCHAR );
    $sth->execute();
  }

  my $exDB;

  while ( my $rowcache = $precache || $sth->fetchall_arrayref( undef, $max_rows ) ) {
    $precache = undef;
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $dbID,                $dbprimaryId,
           $displayid,           $version,
           $priority,
           $dbname,              $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      if ( !defined($exDB) ) {
        $exDB =
          Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $dbname,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
			   -description         => $description
          );


      }

      if ($synonym) { $exDB->add_synonym($synonym) }

    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return $exDB;
} ## end sub fetch_by_db_accession


=head2 store

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbEntry
               The DBEntry (xref) to be stored
  Arg [2]    : Int $ensID
               The dbID of an EnsEMBL object to associate with this external
               database entry
  Arg [3]    : string $ensType ('Transcript', 'Translation', 'Gene')
               The type of EnsEMBL object that this external database entry is
               being associated with.
  Arg [4]    : boolean $ignore_release
               If unset or zero, will require that the release string
               of the DBEntry object is identical to the release of the
               external database.  If set and non-zero, will ignore the
               release information.
  Example    : $dbea->store($db_entry, $transcript_id, 'Transcript');
  Description: Stores a reference to an external database (if it is not stored
               already) and associates an EnsEMBL object of a specified type
               with the external identifier.
  Returntype : int - the dbID of the newly created external refernce
  Exceptions : thrown when invalid dbID is passed to this method
  Caller     : scripts which load Xrefs and ObjectXrefs, etc. into Ensembl
  Status     : Stable

=cut

sub store {
  my ( $self, $dbEntry, $ensID, $ensType, $ignore_release, $master_xref ) = @_;

  my $dbJustInserted;

  #
  # backwards compatibility check:
  # check if $ensID is an object; if so, use $obj->dbID
  #
  my $ensembl_id;

  if ( defined($ensID) ) {
    if ( $ensID =~ /^\d+$/ ) {
      $ensembl_id = $ensID;
    } elsif (    ref($ensID) eq 'Bio::EnsEMBL::Gene'
              or ref($ensID) eq 'Bio::EnsEMBL::Transcript'
              or ref($ensID) eq 'Bio::EnsEMBL::Translation' 
              or ref($ensID) eq 'Bio::EnsEMBL::OperonTranscript'
              or ref($ensID) eq 'Bio::EnsEMBL::Operon' 
              )
    {
      warning(   "You should pass DBEntryAdaptor->store() "
               . "a dbID rather than an ensembl object "
               . "to store the xref on" );

      if ( defined( $ensID->dbID() ) ) {
        $ensembl_id = $ensID->dbID();
      } else {
        throw( sprintf( "%s %s doesn't have a dbID, can't store xref",
                        $ensType, $ensID->display_id() ) );
      }
    } else {
      throw("Invalid dbID passed to DBEntryAdaptor->store()");
    }
  }
  
  
  
    # Ensure external_db contains a record of the intended xref source
    my $dbRef;
    $dbRef = $self->_check_external_db($dbEntry,$ignore_release);

    # Attempt to insert DBEntry
    my $xref_id = $self->_store_or_fetch_xref($dbEntry,$dbRef);
    $dbEntry->dbID($xref_id); #keeps DBEntry in sync with database
    ### Attempt to create an object->xref mapping
    my $object_xref_id;
    if ($ensembl_id) { $object_xref_id = $self->_store_object_xref_mapping($ensembl_id,$dbEntry,$ensType, $ignore_release)};

    if (defined $master_xref && defined $object_xref_id) { $self->_store_dependent_xref_mapping($object_xref_id, $dbEntry, $master_xref); }
    
    return $xref_id;
}




=head2 update

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbentry
               The dbentry to update
  Example    : $dbentry_adaptor->update($dbentry);
  Description: Updates the dbprimary_acc, display_label, version, description, info_type
               and info_text of a dbentry in the database.
  Returntype : None
  Exceptions : thrown if the $dbentry is not a Bio::EnsEMBL::DBEntry
  Caller     : general
  Status     : Stable

=cut

sub update {
  my ($self, $dbEntry) = @_;

  if (!defined $dbEntry || !ref $dbEntry || !$dbEntry->isa('Bio::EnsEMBL::DBEntry')) {
    throw("Must update a dbentry object, not a $dbEntry");
  }

  my $update_dbentry_sql = qq(
     UPDATE xref
        SET dbprimary_acc = ?,
            display_label = ?,
            version = ?,
            description = ?,
            external_db_id = ?,
            info_type = ?,
            info_text = ?
      WHERE xref_id = ?
  );

  my $dbRef = $self->_check_external_db($dbEntry);
  my $display_id = $dbEntry->display_id;
  $display_id = '' unless defined $display_id; # SQLite doesn't ignore NOT NULL errors

  my $sth = $self->prepare($update_dbentry_sql);

  $sth->bind_param(1, $dbEntry->primary_id,SQL_VARCHAR);
  $sth->bind_param(2, $display_id,SQL_VARCHAR);
  $sth->bind_param(3, ($dbEntry->version || q{0}),SQL_VARCHAR);
  $sth->bind_param(4, $dbEntry->description,SQL_VARCHAR);
  $sth->bind_param(5, $dbRef,SQL_INTEGER);
  $sth->bind_param(6, ($dbEntry->info_type || 'NONE'), SQL_VARCHAR);
  $sth->bind_param(7, ($dbEntry->info_text || ''), SQL_VARCHAR);

  $sth->bind_param(8, $dbEntry->dbID(), SQL_INTEGER);

  $sth->execute();

} ## end sub update
    
sub _store_object_xref_mapping {
    my $self = shift;
    my $ensembl_id = shift;
    my $dbEntry = shift;
    my $ensembl_type = shift;
    my $ignore_release = shift;
    
    my $dbc = $self->dbc();
    my $sql_helper = $dbc->sql_helper();
    
    if (not defined ($ensembl_type)) { warning("No Ensembl data type provided for new xref");}
    
    my $analysis_id;
    if ( $dbEntry->analysis() ) {
        $analysis_id = $self->db()->get_AnalysisAdaptor->store( $dbEntry->analysis() );
    } else {
        $analysis_id = 0; ## This used to be undef, but uniqueness in mysql requires a value
    }
    
    my $insert_ignore = $self->insert_ignore_clause();
    my $insert_params = [
      [$dbEntry->dbID(),              SQL_INTEGER],
      [$ensembl_type,                 SQL_VARCHAR],
      [$ensembl_id,                   SQL_INTEGER],
      [$dbEntry->linkage_annotation(),SQL_VARCHAR],
      [$analysis_id,                  SQL_INTEGER],
    ];
    my $base_object_xref_insert_sql = qq(
      ${insert_ignore} INTO object_xref
            ( xref_id,
              ensembl_object_type,
              ensembl_id,
              linkage_annotation,
              analysis_id )
        VALUES ( ?, ?, ?, ?, ? )
    );
    
    my $object_xref_id;
    
    #If MySQL we can optimise into a single insert query to deal with object_xref_id retrieval on insert ignore
    if($dbc->driver() eq 'mysql') {
      my $sql = $base_object_xref_insert_sql.'ON DUPLICATE KEY UPDATE object_xref_id=LAST_INSERT_ID(object_xref_id)';
      $sql_helper->execute_update(-SQL => $sql, -PARAMS => $insert_params);
      $object_xref_id = $self->last_insert_id('object_xref_id', undef, 'object_xref');
    }
    #Otherwise we attempt & retrieve on failure
    else {
      my $updated_rows = $sql_helper->execute_update(-SQL => $base_object_xref_insert_sql, -PARAMS => $insert_params);
      if($updated_rows == 1) {
        $object_xref_id = $self->last_insert_id('object_xref_id', undef, 'object_xref');
      }
      else {
        my $sql = 'select object_xref_id from object_xref where xref_id =? and ensembl_object_type =? and ensembl_id =? and analysis_id =?';
        my $params = [
          [$dbEntry->dbID(),  SQL_INTEGER],
          [$ensembl_type,     SQL_VARCHAR],
          [$ensembl_id,       SQL_INTEGER],
          [$analysis_id,      SQL_INTEGER],
        ];
        $object_xref_id = $sql_helper->execute_single_result(-SQL => $sql, -PARAMS => $params);
      }
    }
    
    $dbEntry->adaptor($self); # hand Adaptor to dbEntry for future use with OntologyXrefs
    # This is here because everything else assumes it can reuse the $sth scalar ... grrr
    my $sth;
        if ( $dbEntry->isa('Bio::EnsEMBL::IdentityXref') ) {
        # Can attempt multiple inserts as object_xref_id is the PRIMARY KEY
        $sth = $self->prepare( "
             ${insert_ignore} INTO identity_xref
           ( object_xref_id,
             xref_identity,
             ensembl_identity,
             xref_start,
             xref_end  ,
             ensembl_start,
             ensembl_end,
             cigar_line,
             score,
             evalue )
           VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )" );
        $sth->bind_param( 1, $object_xref_id,            SQL_INTEGER );
        $sth->bind_param( 2, $dbEntry->xref_identity,    SQL_INTEGER );
        $sth->bind_param( 3, $dbEntry->ensembl_identity, SQL_INTEGER );
        $sth->bind_param( 4, $dbEntry->xref_start,       SQL_INTEGER );
        $sth->bind_param( 5, $dbEntry->xref_end,         SQL_INTEGER );
        $sth->bind_param( 6, $dbEntry->ensembl_start,    SQL_INTEGER );
        $sth->bind_param( 7, $dbEntry->ensembl_end,      SQL_INTEGER );
        $sth->bind_param( 8,  $dbEntry->cigar_line,  SQL_LONGVARCHAR );
        $sth->bind_param( 9,  $dbEntry->score,            SQL_DOUBLE );
        $sth->bind_param( 10, $dbEntry->evalue,           SQL_DOUBLE );
        $sth->execute();
      } elsif ( $dbEntry->isa('Bio::EnsEMBL::OntologyXref') ) {
        $sth = $self->prepare( "
             ${insert_ignore} INTO ontology_xref
                  ( object_xref_id,
                    source_xref_id,
                    linkage_type    )
               VALUES ( ?, ?, ? )" );
        foreach my $info ( @{ $dbEntry->get_all_linkage_info() } ) {
            my ( $linkage_type, $sourceXref ) = @{$info};
            my $sourceXid = undef;
            if ($sourceXref) {
              $sourceXref->is_stored( $self->dbc ) || $self->store($sourceXref, undef, undef, $ignore_release);
              $sourceXid = $sourceXref->dbID;
            }
            $sth->bind_param( 1, $object_xref_id, SQL_INTEGER );
            $sth->bind_param( 2, $sourceXid, SQL_INTEGER );
            $sth->bind_param( 3, $linkage_type,  SQL_VARCHAR );
            $sth->execute();
        } #end foreach
        
        $sth = $self->prepare( "
             ${insert_ignore} INTO associated_xref
                  ( object_xref_id,
                    xref_id,
                    source_xref_id,
                    condition_type,
                    associated_group_id,
                    rank                 )
           VALUES ( ?, ?, ?, ?, ?, ? ) " );
        
        my $annotext = $dbEntry->get_all_associated_xrefs();
        foreach my $ax_group (sort keys %{ $annotext }) {
          my $group = $annotext->{$ax_group};
          my $gsth = $self->prepare( " 
                  INSERT INTO associated_group 
                    ( description )
                  VALUES ( ? )" );
          $gsth->bind_param( 1, $ax_group,     SQL_VARCHAR );
          $gsth->execute();
          my $associatedGid = $self->last_insert_id('associated_group_id', undef, 'associated_group');
          
          foreach my $ax_rank (sort keys %{ $group }) {
            my @ax = @{ $group->{$ax_rank} };
            
            my $associatedXid = undef;
            my $sourceXid = undef;
            
            if (!$ax[0]->dbID) {
              $self->store($ax[0]);
            }
            $associatedXid = $ax[0]->dbID;

            if (!$ax[1]->dbID) {
              $self->store($ax[1]);
            }
            $sourceXid = $ax[1]->dbID;
            
            if (!defined $associatedXid || !defined $sourceXid) {
              next;
            }
            $sth->bind_param( 1, $object_xref_id,     SQL_INTEGER );
            $sth->bind_param( 2, $associatedXid,      SQL_INTEGER );
            $sth->bind_param( 3, $sourceXid,          SQL_INTEGER );
            $sth->bind_param( 4, $ax[2],              SQL_VARCHAR );
            $sth->bind_param( 5, $associatedGid,      SQL_VARCHAR );
            $sth->bind_param( 6, $ax_rank,            SQL_INTEGER );
            $sth->execute();
          }
        } #end foreach
      } #end elsif
    return $object_xref_id;
}


sub _store_dependent_xref_mapping {
  my $self = shift;
  my $object_xref_id = shift;
  my $dbEntry = shift;
  my $master_xref = shift;

  my $insert_ignore = $self->insert_ignore_clause();

  my $sth = $self->prepare( "
    ${insert_ignore} INTO dependent_xref
       ( object_xref_id,
       master_xref_id,
       dependent_xref_id )
       VALUES ( ?, ?, ?)" );
  $sth->bind_param( 1, $object_xref_id,       SQL_INTEGER );
  $sth->bind_param( 2, $master_xref->dbID,    SQL_INTEGER );
  $sth->bind_param( 3, $dbEntry->dbID,        SQL_INTEGER );

  $sth->execute();

}

=head2 get_external_db_id

  Arg [1]    : String
               The external DB name to query by. Supports LIKE statements
  Arg [2]    : String (optional) 
               External DB release to use. If not specified then we 
               will search for NULL db_release entries
  Arg [3]    : Boolean (optional)
               If true we will never look at the db_release value 
               when querying for an external db id
  Description: Looks for the internal identifier of an external DB. You can
               search using direct equality or using like statements specify. We
               only return one value from this method.

               If you want more than one entry use get_external_db_ids().
  Exceptions : None
  Returntype : Int 

=cut

sub get_external_db_id {
  my $self = shift;
  my $db_ids = $self->get_external_db_ids(@_);
  return shift @$db_ids;
}

=head2 get_external_db_ids

  Arg [1]    : String
               The external DB name to query by. Supports LIKE statements
  Arg [2]    : String (optional) 
               External DB release to use. If not specified then we 
               will search for NULL db_release entries
  Arg [3]    : Boolean (optional)
               If true we will never look at the db_release value 
               when querying for an external db id
  Description: Looks for the internal identifier of an external DB. You can
               search using direct equality or using like statements specify. We
               only return one value from this method. Returns more than one value
  Exceptions : None
  Returntype : ArrayRef of Int 

=cut

sub get_external_db_ids {
    my ($self, $db_name, $db_release, $ignore_release) = @_;
    my $sql_helper = $self->dbc->sql_helper;
    my $sql = 'SELECT external_db_id FROM external_db WHERE db_name';
    if ($db_name =~ /%|_/) { 
        $sql .= ' LIKE ?';
    } else {
        $sql .= ' = ?';
    }
    my @bound_params;
    push @bound_params,$db_name;
    unless ($ignore_release) {
        if ($db_release) {
            $sql .= ' AND db_release = ?';
            push @bound_params,$db_release;
        } 
        else {
            $sql .= ' AND db_release is NULL';
        }
    }
    
    my @db_ids = @{ $sql_helper->execute_simple(-SQL => $sql, -PARAMS => \@bound_params) };
    return \@db_ids;
}

=head2 get_distinct_external_dbs

  Description: Queries the external_db table for all unique 
               DB names
  Exceptions : None
  Returntype : ArrayRef of String

=cut

sub get_distinct_external_dbs {
  my ($self) = @_;
  my $sql = 'SELECT DISTINCT db_name FROM external_db';
  my $names =  $self->dbc->sql_helper->execute_simple(-SQL => $sql);
  return [ sort @{$names} ];
}

=head2 _check_external_db 

  Arg [1]    : DBEntry object
  Arg [2]    : Ignore version flag
  Description: Looks for a record of the given external database
  Exceptions : Throws on missing external database entry
  Returntype : Int 

=cut

sub _check_external_db {
  my ($self,$db_entry,$ignore) = @_;
  my ($db_name,$db_release);
  $db_name = $db_entry->dbname();
  $db_release = $db_entry->release();
  my $db_id = $self->get_external_db_id($db_name, $db_release, $ignore);
  
  if ($db_id) {
    return $db_id;
  }
  else {
    throw( sprintf( "external_db [%s] release [%s] does not exist",
                   $db_name, $db_release)
    );
  }
}

=head2 _store_or_fetch_xref

    Arg [1]    : DBEntry object
    Arg [2]    : Database accession for external database
    Description: Thread-safe method for adding xrefs, or otherwise returning
                 an xref ID for the inserted or retrieved xref. Also inserts
                 synonyms for that xref when entire new 
    Returns    : Int - the DB ID of the xref after insertion 
=cut

sub _store_or_fetch_xref {
    my $self = shift;
    my $dbEntry = shift;
    my $dbRef = shift;
    my $xref_id;

    my $display_id = $dbEntry->display_id;
    $display_id = '' unless defined $display_id; # SQLite doesn't ignore NOT NULL errors

    my $insert_ignore = $self->insert_ignore_clause();
    my $sth = $self->prepare( "
       ${insert_ignore} INTO xref
         ( dbprimary_acc,
           display_label,
           version,
           description,
           external_db_id,
           info_type,
           info_text )
         VALUES ( ?, ?, ?, ?, ?, ?, ? ) ");
    $sth->bind_param(1, $dbEntry->primary_id,SQL_VARCHAR);
    $sth->bind_param(2, $display_id,SQL_VARCHAR);
    $sth->bind_param(3, ($dbEntry->version || q{0}),SQL_VARCHAR);
    $sth->bind_param(4, $dbEntry->description,SQL_VARCHAR);
    $sth->bind_param(5, $dbRef,SQL_INTEGER);
    $sth->bind_param(6, ($dbEntry->info_type || 'NONE'), SQL_VARCHAR);
    $sth->bind_param(7, ($dbEntry->info_text || ''), SQL_VARCHAR);

    my $count = $sth->execute();
    $xref_id = $self->last_insert_id('xref_id',undef,'xref') if $count > 0;
    $sth->finish();
    
    if ($xref_id) { #insert was successful, store supplementary synonyms
        # thread safety no longer an issue.
        my $synonym_check_sth = $self->prepare(
                  "SELECT xref_id, synonym
                   FROM external_synonym
                   WHERE xref_id = ?
                   AND synonym = ?");
    
        my $synonym_store_sth = $self->prepare(
            "${insert_ignore} INTO external_synonym
             ( xref_id, synonym ) VALUES ( ?, ? ) ");
    
        my $synonyms = $dbEntry->get_all_synonyms();
        foreach my $syn ( @$synonyms ) {
            $synonym_check_sth->bind_param(1,$xref_id,SQL_INTEGER);
            $synonym_check_sth->bind_param(2,$syn,SQL_VARCHAR);
            $synonym_check_sth->execute();
            my ($dbSyn) = $synonym_check_sth->fetchrow_array();
            $synonym_store_sth->bind_param(1,$xref_id,SQL_INTEGER);
            $synonym_store_sth->bind_param(2,$syn,SQL_VARCHAR);
            $synonym_store_sth->execute() if(!$dbSyn);
        }
        $synonym_check_sth->finish();
        $synonym_store_sth->finish();
        
    } else { # xref_id already exists, retrieve it according to fields in the unique key
        my $sql = 'SELECT xref_id FROM xref 
            WHERE dbprimary_acc = ?
            AND version =?
            AND external_db_id = ?
            AND info_type = ?
            AND info_text = ?';
        my $info_type = $dbEntry->info_type() || 'NONE';
        my $info_text = $dbEntry->info_text() || q{};
        my $version = $dbEntry->version() || q{0};
        $sth = $self->prepare( $sql );

        $sth->bind_param(1, $dbEntry->primary_id,SQL_VARCHAR);
        $sth->bind_param(2, $version, SQL_VARCHAR);
        $sth->bind_param(3, $dbRef, SQL_INTEGER);
        $sth->bind_param(4, $info_type, SQL_VARCHAR);
        $sth->bind_param(5, $info_text, SQL_VARCHAR);
        $sth->execute();



        ($xref_id) = $sth->fetchrow_array();
        $sth->finish;
        if(!$xref_id) {
          my $msg = 'Cannot find an xref id for %s (version=%d) with external db id %d.';
          throw(sprintf($msg, $dbEntry->primary_id(), $version, $dbRef))
        }
    }
    
    return $xref_id;
}

=head2 exists

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
  Example    : if($dbID = $db_entry_adaptor->exists($dbe)) { do stuff; }
  Description: Returns the db id of this DBEntry if it exists in this database
               otherwise returns undef.  Exists is defined as an entry with
               the same external_db and display_id
  Returntype : int
  Exceptions : thrown on incorrect args
  Caller     : GeneAdaptor::store, TranscriptAdaptor::store
  Status     : Stable

=cut

sub exists {
  my ($self, $dbe) = @_ ;

  unless($dbe && ref $dbe && $dbe->isa('Bio::EnsEMBL::DBEntry')) {
    throw("arg must be a Bio::EnsEMBL::DBEntry not [$dbe]");
  }

  my $sth = $self->prepare('SELECT x.xref_id
                            FROM   xref x, external_db xdb
                            WHERE  x.external_db_id = xdb.external_db_id
                            AND    x.display_label = ?
                            AND    xdb.db_name = ?
                            AND    x.dbprimary_acc = ?');

  $sth->bind_param(1,$dbe->display_id,SQL_VARCHAR);
  $sth->bind_param(2,$dbe->dbname,SQL_VARCHAR);
  $sth->bind_param(3,$dbe->primary_id,SQL_VARCHAR);
  $sth->execute();

  my ($dbID) = $sth->fetchrow_array;

  $sth->finish;

  return $dbID;
}


=head2 fetch_all_by_Gene

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               (The gene to retrieve DBEntries for)
  Arg [2]    : optional external database name. SQL wildcards are accepted
  Arg [3]    : optional external_db type. SQL wildcards are accepted
  Example    : @db_entries = @{$db_entry_adaptor->fetch_all_by_Gene($gene)};
  Description: This returns a list of DBEntries associated with this gene.
               Note that this method was changed in release 15.  Previously
               it set the DBLinks attribute of the gene passed in to contain
               all of the gene, transcript, and translation xrefs associated
               with this gene.
  Returntype : listref of Bio::EnsEMBL::DBEntries; may be of type IdentityXref if
               there is mapping data, or OntologyXref if there is linkage data.
  Exceptions : thows if gene object not passed
  Caller     : Bio::EnsEMBL::Gene
  Status     : Stable

=cut

sub fetch_all_by_Gene {
  my ( $self, $gene, $ex_db_reg, $exdb_type ) = @_;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Bio::EnsEMBL::Gene argument expected.");
  }

  return $self->_fetch_by_object_type($gene->dbID(), 'Gene', $ex_db_reg, $exdb_type);
}

=head2 fetch_all_by_Operon

  Arg [1]    : Bio::EnsEMBL::Operon $operon
               (The operon to retrieve DBEntries for)
  Arg [2]    : optional external database name. SQL wildcards are accepted
  Arg [3]    : optional external_db type. SQL wildcards are accepted
  Example    : @db_entries = @{$db_entry_adaptor->fetch_all_by_Operon($operon)};
  Description: This returns a list of DBEntries associated with this operon.
  Returntype : listref of Bio::EnsEMBL::DBEntries; may be of type IdentityXref if
               there is mapping data, or OntologyXref if there is linkage data.
  Exceptions : thows if operon object not passed
  Caller     : general

=cut

sub fetch_all_by_Operon {
  my ( $self, $gene, $ex_db_reg, $exdb_type ) = @_;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Operon')) {
    throw("Bio::EnsEMBL::Operon argument expected.");
  }

  return $self->_fetch_by_object_type($gene->dbID(), 'Operon', $ex_db_reg, $exdb_type);
}


=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript
  Arg [2]    : optional external database name. SQL wildcards are accepted
  Arg [3]    : optional external_db type. SQL wildcards are accepted
  Example    : @db_entries = @{$db_entry_adaptor->fetch_all_by_Transcript($trans)};
  Description: This returns a list of DBEntries associated with this
               transcript. Note that this method was changed in release 15.
               Previously it set the DBLinks attribute of the gene passed in
               to contain all of the gene, transcript, and translation xrefs
               associated with this gene.
  Returntype : listref of Bio::EnsEMBL::DBEntries; may be of type IdentityXref if
               there is mapping data, or OntologyXref if there is linkage data.
  Exceptions : throes if transcript argument not passed
  Caller     : Bio::EnsEMBL::Transcript
  Status     : Stable

=cut

sub fetch_all_by_Transcript {
  my ( $self, $trans, $ex_db_reg, $exdb_type ) = @_;

  if(!ref($trans) || !$trans->isa('Bio::EnsEMBL::Transcript')) {
    throw("Bio::EnsEMBL::Transcript argument expected.");
  }

  return $self->_fetch_by_object_type( $trans->dbID(), 'Transcript', $ex_db_reg, $exdb_type);
}


=head2 fetch_all_by_Translation

  Arg [1]    : Bio::EnsEMBL::Translation $trans
               (The translation to fetch database entries for)
  Arg [2]    : optional external database name. SQL wildcards are accepted
  Arg [3]    : optional externaldb type. SQL wildcards are accepted
  Example    : @db_entries = @{$db_entry_adptr->fetch_all_by_Translation($trans)};
  Description: Retrieves external database entries for an EnsEMBL translation
  Returntype : listref of Bio::EnsEMBL::DBEntries; may be of type IdentityXref if
               there is mapping data, or OntologyXref if there is linkage data.
  Exceptions : throws if translation object not passed
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Translation {
  my ( $self, $trans, $ex_db_reg, $exdb_type ) = @_;

  if(!ref($trans) || !$trans->isa('Bio::EnsEMBL::Translation')) {
    throw('Bio::EnsEMBL::Translation argument expected.');
  }
  if( ! $trans->dbID ){
    warning( "Cannot fetch_all_by_Translation without a dbID" );
    return [];
  }

  return $self->_fetch_by_object_type( $trans->dbID(), 'Translation', $ex_db_reg, $exdb_type );
}



=head2 remove_from_object

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe - The external reference which
               is to be disassociated from an ensembl object.
  Arg [2]    : Bio::EnsEMBL::Storable $object - The ensembl object the
               external reference is to be disassociated from
  Arg [3]    : string $object_type - The type of the ensembl object.
               E.g. 'Gene', 'Transcript', 'Translation'
  Example    :
               # remove all dbentries from this translation
               foreach my $dbe (@{$translation->get_all_DBEntries()}) {
                 $dbe_adaptor->remove($dbe, $translation, 'Translation');
               }
  Description: Removes an association between an ensembl object and a
               DBEntry (xref).  This does not remove the actual xref from
               the database, only its linkage to the ensembl object.
  Returntype : none
  Exceptions : Throw on incorrect arguments.
               Warning if object or dbentry is not stored in this database.
  Caller     : TranscriptAdaptor::remove, GeneAdaptor::remove,
               TranslationAdaptor::remove
  Status     : Stable

=cut

sub remove_from_object {
  my $self = shift;
  my $dbe  = shift;
  my $object = shift;
  my $object_type = shift;

  if(!ref($dbe) || !$dbe->isa('Bio::EnsEMBL::DBEntry')) {
    throw("Bio::EnsEMBL::DBEntry argument expected.");
  }

  if(!ref($object) || !$dbe->isa('Bio::EnsEMBL::Storable')) {
    throw("Bio::EnsEMBL::Storable argument expected.");
  }

  if(!$object_type) {
    throw("object_type string argument expected.");
  }

  # make sure both the dbentry and the object it is allegedly linked to
  # are stored in this database

  if(!$object->is_stored($self->db())) {
    warning("Cannot remove DBEntries for $object_type " . $object->dbID() .
            ". Object is not stored in this database.");
    return;
  }

  if(!$dbe->is_stored($self->db())) {
    warning("Cannot remove DBEntry ".$dbe->dbID() . ". Is not stored " .
            "in this database.");
    return;
  }

  # obtain the identifier of the link from the object_xref table
  #No need to compare linkage_annotation here
  my $sth = $self->prepare
    ("SELECT ox.object_xref_id " .
     "FROM   object_xref ox ".
     "WHERE  ox.xref_id = ? " .
     "AND    ox.ensembl_id = ? " .
     "AND    ox.ensembl_object_type = ?");
  $sth->bind_param(1,$dbe->dbID,SQL_INTEGER);
  $sth->bind_param(2,$object->dbID,SQL_INTEGER);
  $sth->bind_param(3,$object_type,SQL_VARCHAR);
  $sth->execute();

  my ($ox_id) = $sth->fetchrow_array();

  if(!$sth->rows() == 1) {
    $sth->finish();
    return;
  }

  $sth->finish();

  # delete from the tables which contain additional linkage information

  $sth = $self->prepare("DELETE FROM ontology_xref WHERE object_xref_id = ?");
  $sth->bind_param(1,$ox_id,SQL_INTEGER);
  $sth->execute();
  $sth->finish();
  
  $sth = $self->prepare("DELETE FROM associated_xref WHERE object_xref_id = ?");
  $sth->bind_param(1,$ox_id,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  $sth = $self->prepare("DELETE FROM identity_xref WHERE object_xref_id = ?");
  $sth->bind_param(1,$ox_id,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # delete the actual linkage itself
  $sth = $self->prepare("DELETE FROM object_xref WHERE object_xref_id = ?");
  $sth->bind_param(1,$ox_id,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  return;
}


=head2 _fetch_by_object_type

  Arg [1]    : string $ensID
  Arg [2]    : string $ensType (object type to be returned)
  Arg [3]    : optional $exdbname (external database name)
               (may be an SQL pattern containing '%' which matches any
               number of characters)
  Arg [4]    : optional $exdb_type (external database type)
               (may be an SQL pattern containing '%' which matches any
               number of characters)
  Example    : $self->_fetch_by_object_type( $translation_id, 'Translation' )
  Description: Fetches DBEntry by Object type
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.


  Returntype : arrayref of DBEntry objects; may be of type IdentityXref if
               there is mapping data, or OntologyXref if there is linkage data.
  Exceptions : none
  Caller     : fetch_all_by_Gene
               fetch_all_by_Translation
               fetch_all_by_Transcript
  Status     : Stable

=cut

sub _fetch_by_object_type {
  my ( $self, $ensID, $ensType, $exdbname, $exdb_type ) = @_;

  my @out;

  if ( !defined($ensID) ) {
    throw("Can't fetch_by_EnsObject_type without an object");
  }

  if ( !defined($ensType) ) {
    throw("Can't fetch_by_EnsObject_type without a type");
  }

  ### TODO - SCHEMA VERSION HACK. Please remove sometime in the future once we do not want backwards compatibility
  my $schema_version = $self->schema_version();

  #  my $sth = $self->prepare("
  my $columns_sql = (<<COLUMNS_SQL);
    SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label, xref.version,
           exDB.priority,
           exDB.db_name, exDB.db_release, exDB.status, exDB.db_display_name,
           exDB.secondary_db_name, exDB.secondary_db_table,
           oxr.object_xref_id,
           es.synonym,
           idt.xref_identity, idt.ensembl_identity, idt.xref_start,
           idt.xref_end, idt.ensembl_start, idt.ensembl_end,
           idt.cigar_line, idt.score, idt.evalue, oxr.analysis_id,
           gx.linkage_type,
           xref.info_type, xref.info_text, exDB.type, gx.source_xref_id,
           oxr.linkage_annotation, xref.description
COLUMNS_SQL

  my $tables_sql = <<'TABLES_SQL';
    FROM   (xref xref, external_db exDB, object_xref oxr)
    LEFT JOIN external_synonym es on es.xref_id = xref.xref_id
    LEFT JOIN identity_xref idt on idt.object_xref_id = oxr.object_xref_id
    LEFT JOIN ontology_xref gx on gx.object_xref_id = oxr.object_xref_id
TABLES_SQL
  
  my $where_sql = <<'WHERE_SQL';
    WHERE  xref.xref_id = oxr.xref_id
      AND  xref.external_db_id = exDB.external_db_id
      AND  oxr.ensembl_id = ?
      AND  oxr.ensembl_object_type = ?
WHERE_SQL

  if($schema_version >= 72) {
    $columns_sql .= ','; #need the extra comma
    $columns_sql .= <<'EXTRA';
    ax.xref_id, ax.source_xref_id, ax.condition_type, ax.associated_group_id, ax.rank
EXTRA
    $tables_sql .= <<'EXTRA';
    LEFT JOIN associated_xref ax ON ax.object_xref_id = oxr.object_xref_id
    LEFT JOIN associated_group ag ON ax.associated_group_id = ag.associated_group_id
EXTRA

  }

  #Join the above together into one statement
  my $sql = "$columns_sql
$tables_sql
$where_sql";
  # warn $sql; # uncomment me to see the full SQL generated

  if ( defined($exdbname) ) {
    if ( index( $exdbname, '%' ) != -1 ) {
      $sql .= " AND exDB.db_name LIKE "
        . $self->dbc()->db_handle()->quote( $exdbname, SQL_VARCHAR );
    } else {
      $sql .= " AND exDB.db_name = "
        . $self->dbc()->db_handle()->quote( $exdbname, SQL_VARCHAR );
    }
  }

  if ( defined($exdb_type) ) {
    if ( index( $exdb_type, '%' ) != -1 ) {
      $sql .= " AND exDB.type LIKE "
        . $self->dbc()->db_handle()->quote( $exdb_type, SQL_VARCHAR );
    } else {
      $sql .= " AND exDB.type = "
        . $self->dbc()->db_handle()->quote( $exdb_type, SQL_VARCHAR );
    }
  }

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $ensID,   SQL_INTEGER );
  $sth->bind_param( 2, $ensType, SQL_VARCHAR );
  $sth->execute();

  my ( %seen, %linkage_types, %synonyms, %associated_xrefs );

  my $max_rows = 1000;

  while ( my $rowcache = $sth->fetchall_arrayref( undef, $max_rows ) ) {
    while ( my $arrRef = shift( @{$rowcache} ) ) {
      my ( $refID,                     $dbprimaryId,
           $displayid,                 $version,
           $priority,
           $dbname,                    $release,
           $exDB_status,               $exDB_db_display_name,
           $exDB_secondary_db_name,    $exDB_secondary_db_table,
           $objid,                     $synonym,
           $xrefid,                    $ensemblid,
           $xref_start,                $xref_end,
           $ensembl_start,             $ensembl_end,
           $cigar_line,                $score,
           $evalue,                    $analysis_id,
           $linkage_type,              $info_type,
           $info_text,                 $type,
           $source_xref_id,            $link_annotation,
           $description,               $associated_xref_id,
           $source_associated_xref_id, $condition_type,
           $associate_group_id,        $associate_group_rank
      ) = @$arrRef;

      my $linkage_key =
        ( $linkage_type || '' ) . ( $source_xref_id || '' );
      
      my $associated_key =
        ( $condition_type || '' )
        . ( $source_associated_xref_id || '' )
        . ( $associated_xref_id || '' );


      my $analysis = undef;
      if ( defined($analysis_id) ) {
	$analysis =
	  $self->db()->get_AnalysisAdaptor()->fetch_by_dbID($analysis_id);
      }

      my %obj_hash = ( 'adaptor'            => $self,
                       'dbID'               => $refID,
                       'primary_id'         => $dbprimaryId,
                       'display_id'         => $displayid,
                       'version'            => $version,
                       'release'            => $release,
                       'info_type'          => $info_type,
                       'info_text'          => $info_text,
                       'type'               => $type,
                       'secondary_db_name'  => $exDB_secondary_db_name,
                       'secondary_db_table' => $exDB_secondary_db_table,
                       'dbname'             => $dbname,
                       'description'        => $description,
                       'linkage_annotation' => $link_annotation,
                       'analysis'           => $analysis,
		       'ensembl_object_type' => $ensType,
		       'ensembl_id'          => $ensID );

      # Using an outer join on the synonyms as well as on identity_xref,
      # we now have to filter out the duplicates (see v.1.18 for
      # original). Since there is at most one identity_xref row per
      # xref, this is easy enough; all the 'extra' bits are synonyms.
      my $source_xref;
      my $associated_xref;
      my $source_associated_xref;
      if ( !$seen{$refID} ) {
	
	my $exDB;
        if ( ( defined($xrefid) ) ) {  # an xref with similarity scores
          $exDB = Bio::EnsEMBL::IdentityXref->new_fast( \%obj_hash );
          $exDB->xref_identity($xrefid);
          $exDB->ensembl_identity($ensemblid);

          $exDB->cigar_line($cigar_line);
          $exDB->xref_start($xref_start);
          $exDB->xref_end($xref_end); # was not here before 14th Jan 2009 ????
          $exDB->ensembl_start($ensembl_start);
          $exDB->ensembl_end($ensembl_end);
          $exDB->score($score);
          $exDB->evalue($evalue);

        } elsif ( defined $linkage_type && $linkage_type ne "" ) {
          $exDB = Bio::EnsEMBL::OntologyXref->new_fast( \%obj_hash );
          $source_xref = ( defined($source_xref_id)
                              ? $self->fetch_by_dbID($source_xref_id)
                              : undef );
          $exDB->add_linkage_type( $linkage_type, $source_xref || () );
          $linkage_types{$refID}->{$linkage_key} = 1;
          
          # Add associated Xref annotations to the OntologyXref entry.
          if ( defined $associated_xref_id && $associated_xref_id ne "" ) {
            $exDB = Bio::EnsEMBL::OntologyXref->new_fast( \%obj_hash );
            $associated_xref = ( defined($associated_xref_id)
                                ? $self->fetch_by_dbID($associated_xref_id)
                                : undef );
            $source_associated_xref = ( defined($source_associated_xref_id)
                                ? $self->fetch_by_dbID($source_associated_xref_id)
                                : undef );
            if ( defined($associated_xref) ) {
              my $ct = '';
              if ( defined $condition_type ) {
                $ct = $condition_type;
              }
              $exDB->add_linked_associated_xref( $associated_xref, $source_associated_xref, $ct, $associate_group_id, $associate_group_rank );
            }
          }

        } else {
          $exDB = Bio::EnsEMBL::DBEntry->new_fast( \%obj_hash );
        }

        if ( defined($exDB_status) ) { $exDB->status($exDB_status) }

        $exDB->priority($priority);
        $exDB->db_display_name($exDB_db_display_name);

        push( @out, $exDB );
        $seen{$refID} = $exDB;

      } ## end if ( !$seen{$refID} )

      # $exDB still points to the same xref, so we can keep adding GO
      # evidence tags or synonyms.

      if ( defined($synonym) && !$synonyms{$refID}->{$synonym} ) {
        if ( defined($synonym) ) {
          $seen{$refID}->add_synonym($synonym);
        }
        $synonyms{$refID}->{$synonym} = 1;
      }

      if (    defined($linkage_type)
           && $linkage_type ne ""
           && !$linkage_types{$refID}->{$linkage_key} )
      {
        $source_xref = ( defined($source_xref_id)
                            ? $self->fetch_by_dbID($source_xref_id)
                            : undef );
        $seen{$refID}
          ->add_linkage_type( $linkage_type, $source_xref || () );
        $linkage_types{$refID}->{$linkage_key} = 1;
      }
      
      if (    defined($associated_xref_id)
           && $associated_xref_id ne ""
           && !$associated_xrefs{$refID}->{$associated_key} )
      {
        $associated_xref = ( defined($associated_xref_id)
                            ? $self->fetch_by_dbID($associated_xref_id)
                            : undef );
        $source_associated_xref = ( defined($source_associated_xref_id)
                                ? $self->fetch_by_dbID($source_associated_xref_id)
                                : undef );
        if ( defined($associated_xref) ) {
          my $ct = '';
          if ( defined $condition_type ) {
            $ct = $condition_type;
          }
          $seen{$refID}->add_linked_associated_xref( $associated_xref, $source_associated_xref, $ct, $associate_group_id, $associate_group_rank );
        }
        
        $linkage_types{$refID}->{$linkage_key} = 1;
      }
      
    } ## end while ( my $arrRef = shift...
  } ## end while ( my $rowcache = $sth...

  return \@out;
} ## end sub _fetch_by_object_type

=head2 list_gene_ids_by_external_db_id

  Arg [1]    : string $external_id
  Example    : @gene_ids = $dbea->list_gene_ids_by_external_db_id(1020);
  Description: Retrieve a list of geneid by an external identifier that
               is linked to any of the genes transcripts, translations
               or the gene itself.
               NOTE:  If more than one external identifier has the
               same primary accession then genes for each of these is
               returned.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_gene_ids_by_external_db_id {
   my ($self,$external_db_id, $linkage_type) = @_;

   my %T = map { ($_, 1) }
       $self->_type_by_external_db_id( $external_db_id, 'Translation', 'gene', $linkage_type ),
       $self->_type_by_external_db_id( $external_db_id, 'Transcript',  'gene', $linkage_type ),
       $self->_type_by_external_db_id( $external_db_id, 'Gene', undef, $linkage_type );
   return keys %T;
}

=head2 list_gene_ids_by_extids

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Arg [3]    : Boolean override, see _type_by_external_id
  Example    : @gene_ids = $dbea->list_gene_ids_by_extids('CDPX');
  Description: Retrieve a list of geneid by an external identifier that is
               linked to  any of the genes transcripts, translations or the
               gene itself
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_gene_ids_by_extids {
  my ( $self, $external_name, $external_db_name, $override ) = @_;

  my %T = map { ( $_, 1 ) }
    $self->_type_by_external_id( $external_name, 'Translation', 'gene',
                                 $external_db_name, $override ),
    $self->_type_by_external_id( $external_name, 'Transcript', 'gene',
                                 $external_db_name, $override ),
    $self->_type_by_external_id( $external_name, 'Gene', undef,
                                 $external_db_name, $override );

  return keys %T;
}


=head2 list_transcript_ids_by_extids

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Arg [3]    : Boolean override, see _type_by_external_id
  Example    : @tr_ids = $dbea->list_transcript_ids_by_extids('BCRA2');
  Description: Retrieve a list transcript ids by an external identifier that
               is linked to any of the genes transcripts, translations or the
               gene itself
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_transcript_ids_by_extids {
  my ( $self, $external_name, $external_db_name, $override ) = @_;

  my %T = map { ( $_, 1 ) }
    $self->_type_by_external_id( $external_name, 'Translation',
                                 'transcript',   $external_db_name, $override
    ),
    $self->_type_by_external_id( $external_name, 'Transcript', undef,
                                 $external_db_name, $override );

  return keys %T;
}

sub list_transcript_ids_by_external_db_id {
    my ( $self, $external_db_id, $linkage_type) = @_;
    my %T = map { ( $_, 1 ) }
        $self->_type_by_external_db_id( $external_db_id, 'Translation', 'transcript', $linkage_type ),
        $self->_type_by_external_db_id( $external_db_id, 'Transcript',  'transcript', $linkage_type );
    return keys %T;
}

=head2 list_translation_ids_by_extids

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Arg [3]    : Boolean override, see _type_by_external_id
  Example    : @tr_ids = $dbea->list_translation_ids_by_extids('GO:0004835');
  Description: Gets a list of translation IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_translation_ids_by_extids {
  my ( $self, $external_name, $external_db_name, $override ) = @_;

  return
    $self->_type_by_external_id( $external_name, 'Translation', undef,
                                 $external_db_name, $override );
}

sub list_translation_ids_by_external_db_id {
    my ( $self, $external_db_id, $linkage_type) = @_;
    return $self->_type_by_external_db_id( $external_db_id, 'Translation', undef, $linkage_type ),
}

=head2 _type_by_external_id

  Arg [1]    : string $name - dbprimary_acc
  Arg [2]    : string $ensType - ensembl_object_type
  Arg [3]    : (optional) string $extraType
  Arg [4]    : (optional) string $external_db_name
  	           other object type to be returned
  Arg [5]    : Boolean override to force _ to be treated as an SQL 'any'
               This is usually optimised out for query speed due to 
               large numbers of names like NM_00...
  Example    : $self->_type_by_external_id($name, 'Translation');
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
               SQL wildcards can be used in the external id, 
               but overly generic queries (two characters) will be prevented.
  Description: Gets
  Returntype : list of dbIDs (gene_id, transcript_id, etc.)
  Exceptions : none
  Caller     : list_translation_ids_by_extids
               translationids_by_extids
  	           geneids_by_extids
  Status     : Stable

=cut

sub _type_by_external_id {
  my ( $self, $name, $ensType, $extraType, $external_db_name, $override ) = @_;

  # $name has SQL wildcard support
  # = or LIKE put into SQL statement, and open queries like % or A% are rejected.
  my $comparison_operator;
  if ($name =~ /[_%\[]/ ) {
    $comparison_operator = "LIKE";
    if ($name =~ /^.?%/ && !$override) {
      warn "External $ensType name $name is too vague and will monopolise database resources. Please use a more specific $ensType name.\n";
      return;
    }
    elsif ($name =~ /^\w\w_/ && !$override) {
        # For entries such as NM_00000065, escape the _ so that SQL LIKE does not have to scan entire table
        # Escape only the _ in the third character position
        $name =~ s/(?<=\w\w)(?=_)/\\/;
    }
  }
  else {
    $comparison_operator = "=";
  }


  my $from_sql  = '';
  my $where_sql = '';
  my $ID_sql    = 'oxr.ensembl_id';

  if ( defined($extraType) ) {
    if ( lc($extraType) eq 'translation' ) {
      $ID_sql = 'tl.translation_id';
    } else {
      $ID_sql = "t.${extraType}_id";
    }

    if ( lc($ensType) eq 'translation' ) {
      $from_sql  = 'transcript t, translation tl, ';
      $where_sql = qq(
          t.transcript_id = tl.transcript_id AND
          tl.translation_id = oxr.ensembl_id AND
          t.is_current = 1 AND
      );
    } else {
      $from_sql  = 'transcript t, ';
      $where_sql = 't.'
        . lc($ensType)
        . '_id = oxr.ensembl_id AND '
        . 't.is_current = 1 AND ';
    }
  }
  
  my $multispecies = $self->db()->is_multispecies();

  if ( lc($ensType) eq 'gene' ) {
    $from_sql = 'gene g, ';
    $from_sql .= 'seq_region s, coord_system cs, ' if $multispecies;
    
    $where_sql = 'g.gene_id = oxr.ensembl_id AND g.is_current = 1 AND ';
    if($multispecies) {
      $where_sql .= <<'SQL';
g.seq_region_id = s.seq_region_id AND
s.coord_system_id = cs.coord_system_id AND
cs.species_id = ? AND 
SQL
    }
  } 
  elsif ( lc($ensType) eq 'transcript' ) {
    $from_sql = 'transcript t, ';
    $from_sql .= 'seq_region s, coord_system cs, ' if $multispecies;
    
    $where_sql = 't.transcript_id = oxr.ensembl_id AND t.is_current = 1 AND ';
    if($multispecies) {
      $where_sql .= <<'SQL';
t.seq_region_id = s.seq_region_id AND
s.coord_system_id = cs.coord_system_id AND
cs.species_id = ? AND 
SQL
    }
  } 
  elsif ( lc($ensType) eq 'translation' ) {
    $from_sql = 'translation tl, transcript t, ';
    $from_sql .= 'seq_region s, coord_system cs, ' if $multispecies;

    $where_sql = 't.transcript_id = tl.transcript_id AND tl.translation_id = oxr.ensembl_id AND t.is_current = 1 AND ';
    if($multispecies) {
      $where_sql .= <<'SQL';
t.seq_region_id = s.seq_region_id AND
s.coord_system_id = cs.coord_system_id AND
cs.species_id = ? AND 
SQL
    }
  }

  if ( defined($external_db_name) ) {
    # Involve the 'external_db' table to limit the hits to a particular
    # external database.

    $from_sql .= 'external_db xdb, ';
    $where_sql .=
        'xdb.db_name LIKE '
      . $self->dbc()->db_handle()->quote( $external_db_name . '%' )
      . ' AND xdb.external_db_id = x.external_db_id AND';
  }

  my @queries;
  push (@queries, qq(
        SELECT    $ID_sql
        FROM      $from_sql
                  xref x,
                  object_xref oxr
        WHERE     $where_sql
                    x.dbprimary_acc $comparison_operator ?
        AND         x.xref_id = oxr.xref_id
        AND         oxr.ensembl_object_type = ?
    ));

  push (@queries, qq(
        SELECT    $ID_sql
        FROM      $from_sql
                  xref x,
                  object_xref oxr
        WHERE     $where_sql
                    x.display_label $comparison_operator ?
        AND         x.xref_id = oxr.xref_id
        AND         oxr.ensembl_object_type = ?
    ));

  if ( defined($external_db_name) ) {
    # If we are given the name of an external database, we need to join
    # between the 'xref' and the 'object_xref' tables on 'xref_id'.

    push (@queries, qq(
      SELECT    $ID_sql
      FROM      $from_sql
                external_synonym syn,
                object_xref oxr,
                xref x
      WHERE     $where_sql
                syn.synonym $comparison_operator ?
      AND       syn.xref_id = oxr.xref_id
      AND       oxr.ensembl_object_type = ?
      AND       x.xref_id = oxr.xref_id));

  } else {
    # If we weren't given an external database name, we can get away
    # with less joins here.

    push (@queries, qq(
      SELECT    $ID_sql
      FROM      $from_sql
                external_synonym syn,
                object_xref oxr
      WHERE     $where_sql
                syn.synonym $comparison_operator ?
      AND       syn.xref_id = oxr.xref_id
      AND       oxr.ensembl_object_type = ?));

  }

  my %result;
  my $h = $self->dbc()->sql_helper();
  my @params = ([$name, SQL_VARCHAR], [$ensType, SQL_VARCHAR]);
  unshift(@params, [$self->species_id(), SQL_INTEGER] ) if $multispecies;
  foreach my $query (@queries) {
    $h->execute_no_return(-SQL => $query, -PARAMS => \@params, -CALLBACK => sub {
      my ($row) = @_;
      my ($id) = @{$row};
      $result{$id} = 1;
    });
  }

  return keys(%result);

} ## end sub _type_by_external_id

=head2 _type_by_external_db_id

  Arg [1]    : integer $type - external_db_id
  Arg [2]    : string $ensType - ensembl_object_type
  Arg [3]    : (optional) string $extraType
  	       other object type to be returned. This references the _id fields of
  	       the transcript table, and if left unset defaults to the translation_id
  	       from the translation table.
  Arg [4]    : (optional) string $linkage_type
               provides a link to ontology_xref and limits by an ontology linkage type
  Example    : $self->_type_by_external_db_id(1030, 'Translation');
  Description: Gets.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : list of dbIDs (gene_id, transcript_id, etc.)
  Exceptions : none
  Caller     : list_translation_ids_by_extids
               translationids_by_extids
  	           geneids_by_extids
  Status     : Stable

=cut

sub _type_by_external_db_id {
  my ($self, $external_db_id, $ensType, $extraType, $linkage_type) = @_;
  throw "No external DB identifer given" unless defined $external_db_id;
  throw "No Ensembl type given" unless defined $ensType;

  my $from_sql = '';
  my $where_sql = '';
  my $ID_sql = "oxr.ensembl_id";
  my $lcEnsType = lc($ensType);

  if (defined $extraType) {
    if (lc($extraType) eq 'translation') {
      $ID_sql = "tl.translation_id";
    } else {
      $ID_sql = "t.${extraType}_id";
    }

    if ($lcEnsType eq 'translation') {
      $from_sql = 'transcript t, translation tl, ';
      $where_sql = qq(
          t.transcript_id = tl.transcript_id AND
          tl.translation_id = oxr.ensembl_id AND
          t.is_current = 1 AND
      );
    } else {
      $from_sql = 'transcript t, ';
      $where_sql = 't.'.$lcEnsType.'_id = oxr.ensembl_id AND '.
          't.is_current = 1 AND ';
    }
  }

  if ($lcEnsType eq 'gene') {
    $from_sql = 'gene g, ';
    $where_sql = 'g.gene_id = oxr.ensembl_id AND g.is_current = 1 AND ';
  } elsif ($lcEnsType eq 'transcript') {
    $from_sql = 'transcript t, ';
    $where_sql = 't.transcript_id = oxr.ensembl_id AND t.is_current = 1 AND ';
  } elsif ($lcEnsType eq 'translation') {
    $from_sql = 'transcript t, translation tl, ';
    $where_sql = qq(
        t.transcript_id = tl.transcript_id AND
        tl.translation_id = oxr.ensembl_id AND
        t.is_current = 1 AND 
    );
  }

  if($linkage_type) {
    $from_sql .= 'ontology_xref ontx, ';
    $where_sql .= 'oxr.object_xref_id = ontx.object_xref_id AND ';
    $where_sql .= 'ontx.linkage_type =? AND';
  }

  my $query =
    "SELECT $ID_sql
       FROM $from_sql xref x, object_xref oxr
      WHERE $where_sql x.external_db_id = ? AND
  	     x.xref_id = oxr.xref_id AND oxr.ensembl_object_type= ?";

  my @params;
  push(@params, [$linkage_type, SQL_VARCHAR]) if $linkage_type;
  push(@params,
    [$external_db_id, SQL_INTEGER],
    [$ensType, SQL_VARCHAR],
  );

  my %result;
  $self->dbc()->sql_helper()->execute_no_return(-SQL => $query, -PARAMS => \@params, -CALLBACK => sub {
     my ($row) = @_;
     my ($id) = @{$row};
    $result{$id} = 1;
  });

  return keys(%result);
} ## end _type_by_external_db_id


=head2 fetch_all_by_description

  Arg [1]    : string description to search for. Include % etc in this string
  Arg [2]    : <optional> string $dbname. Name of the database to search

  Example    : @canc_refs = @{$db_entry_adaptor->fetch_all_by_description("%cancer%")};
               @db_entries = @{$db_entry_adaptor->fetch_all_by_description("%cancer%","MIM_MORBID")};
  Description: Retrieves DBEntries that match the description.
               Optionally you can search on external databases type.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : ref to array of Bio::EnsEMBL::DBSQL::DBEntry
  Exceptions : None.
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_description {
  my ( $self, $description, $dbname ) = @_;

  my @results = ();

  my $sql =
    "SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label,
           xref.version,
           exDB.priority,
           exDB.db_name, exDB.db_display_name, exDB.db_release, es.synonym,
           xref.info_type, xref.info_text, exDB.type, exDB.secondary_db_name,
           exDB.secondary_db_table, xref.description
    FROM   (xref, external_db exDB)
    LEFT JOIN external_synonym es on es.xref_id = xref.xref_id
    WHERE  xref.description like ?
    AND    xref.external_db_id = exDB.external_db_id";

  if ( defined($dbname) ) { $sql .= " AND exDB.db_name = ? " }

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $description, SQL_VARCHAR );

  if ( defined($dbname) ) {
    $sth->bind_param( 2, $dbname, SQL_VARCHAR );
  }

  $sth->execute();

  my $max_rows = 1000;

  while ( my $rowcache = $sth->fetchall_arrayref( undef, $max_rows ) ) {
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $dbID,                $dbprimaryId,
           $displayid,           $version,
           $priority,
           $ex_dbname,           $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      my $exDB =
        Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $ex_dbname,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
                           -description         => $description
        );

      if ($synonym) { $exDB->add_synonym($synonym) }

      push @results, $exDB;

    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return \@results;
} ## end sub fetch_all_by_description


=head2 fetch_all_by_source

  Arg [1]    : string source to search for. Include % etc in this string
               if you want to use SQL patterns

  Example    : @unigene_refs = @{$db_entry_adaptor->fetch_all_by_source("%unigene%")};
  Description: Retrieves DBEntrys that match the source name.
  Returntype : ref to array of Bio::EnsEMBL::DBSQL::DBEntry
  Exceptions : None.
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_source {
  my ( $self, $source ) = @_;

  my @results = ();

  my $sql =
    "SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label,
           xref.version,
           exDB.priority,
           exDB.db_name, exDB.db_display_name, exDB.db_release, es.synonym,
           xref.info_type, xref.info_text, exDB.type, exDB.secondary_db_name,
           exDB.secondary_db_table, xref.description
    FROM   (xref, external_db exDB)
    LEFT JOIN external_synonym es on es.xref_id = xref.xref_id
    WHERE  exDB.db_name like ?
    AND    xref.external_db_id = exDB.external_db_id";


  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $source, SQL_VARCHAR );

  $sth->execute();

  my $max_rows = 1000;

  while ( my $rowcache = $sth->fetchall_arrayref( undef, $max_rows ) ) {
    while ( my $arrayref = shift( @{$rowcache} ) ) {
      my ( $dbID,                $dbprimaryId,
           $displayid,           $version,
           $priority,
           $dbname,              $db_display_name,
           $release,             $synonym,
           $info_type,           $info_text,
           $type,                $secondary_db_name,
           $secondary_db_table,  $description
      ) = @$arrayref;

      my $exDB =
        Bio::EnsEMBL::DBEntry->new(
                           -adaptor             => $self,
                           -dbID                => $dbID,
                           -primary_id          => $dbprimaryId,
                           -display_id          => $displayid,
                           -version             => $version,
                           -release             => $release,
                           -dbname              => $dbname,
                           -priority            => $priority,
                           -db_display_name     => $db_display_name,
                           -info_type           => $info_type,
                           -info_text           => $info_text,
                           -type                => $type,
                           -secondary_db_name   => $secondary_db_name,
                           -secondary_db_table  => $secondary_db_table,
                           -description         => $description
        );

      if ($synonym) { $exDB->add_synonym($synonym) }

      push @results, $exDB;

    } ## end while ( my $arrayref = shift...
  } ## end while ( my $rowcache = $sth...

  $sth->finish();

  return \@results;
} ## end sub fetch_all_by_source


=head2 fetch_all_synonyms

  Arg [1]    : dbID of DBEntry to fetch synonyms for. Used in lazy loading of synonyms.

  Example    : @canc_refs = @{$db_entry_adaptor->fetch_all_synonyms(1234)};
  Description: Fetches the synonyms for a particular DBEntry.
  Returntype : listref of synonyms. List referred to may be empty if there are no synonyms.
  Exceptions : None.
  Caller     : General
  Status     : At Risk

=cut


sub fetch_all_synonyms {
  my ( $self, $dbID ) = @_;

  my @synonyms = ();

  my $sth =
    $self->prepare( "SELECT synonym "
      . "FROM external_synonym "
      . "WHERE xref_id = ?" );

  $sth->bind_param( 1, $dbID, SQL_INTEGER );

  $sth->execute();

  my $synonym;
  $sth->bind_col(1, \$synonym);

  while ( $sth->fetch() ) {
    push( @synonyms, $synonym );
  }

  return \@synonyms;
}


=head2 get_db_name_from_external_db_id

  Arg [1]    : external_dbid of database to get the database_name
  Example    : my $db_name = $db_entry_adaptor->get_db_name_from_external_db_id(1100);
  Description: Gets the database name for a certain external_db_id
  Returntype : scalar
  Exceptions : None.
  Caller     : General
  Status     : At Risk

=cut

sub get_db_name_from_external_db_id{
    my $self = shift;
    my $external_db_id = shift;

    my $sth = $self->prepare("SELECT db_name FROM external_db WHERE external_db_id = ?");

    $sth->bind_param(1, $external_db_id, SQL_INTEGER);
    $sth->execute();
    my ($db_name) = $sth->fetchrow_array();
    $sth->finish();
    return $db_name;

}

1;


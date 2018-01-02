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

Bio::EnsEMBL::ArchiveStableIdAdaptor

=head1 SYNOPSIS

  my $registry = "Bio::EnsEMBL::Registry";

  my $archiveStableIdAdaptor =
    $registry->get_adaptor( 'Human', 'Core', 'ArchiveStableId' );

  my $stable_id = 'ENSG00000068990';

  my $arch_id = $archiveStableIdAdaptor->fetch_by_stable_id($stable_id);

  print("Latest incarnation of this stable ID:\n");
  printf( "  Stable ID: %s.%d\n",
    $arch_id->stable_id(), $arch_id->version() );
  print("  Release: "
      . $arch_id->release() . " ("
      . $arch_id->assembly() . ", "
      . $arch_id->db_name()
      . ")\n" );

  print "\nStable ID history:\n\n";

  my $history =
    $archiveStableIdAdaptor->fetch_history_tree_by_stable_id(
    $stable_id);

  foreach my $a ( @{ $history->get_all_ArchiveStableIds } ) {
    printf( "  Stable ID: %s.%d\n", $a->stable_id(), $a->version() );
    print("  Release: "
        . $a->release() . " ("
        . $a->assembly() . ", "
        . $a->db_name()
        . ")\n\n" );
  }

=head1 DESCRIPTION

ArchiveStableIdAdaptor does all SQL to create ArchiveStableIds and works
of

  stable_id_event
  mapping_session
  peptite_archive
  gene_archive

tables inside the core database.

This whole module has a status of At Risk as it is under development.

=head1 METHODS

  fetch_by_stable_id
  fetch_by_stable_id_version
  fetch_by_stable_id_dbname
  fetch_all_by_archive_id
  fetch_predecessors_by_archive_id
  fetch_successors_by_archive_id
  fetch_history_tree_by_stable_id
  add_all_current_to_history
  list_dbnames
  previous_dbname
  next_dbname
  get_peptide
  get_current_release
  get_current_assembly

=head1 RELATED MODULES

  Bio::EnsEMBL::ArchiveStableId
  Bio::EnsEMBL::StableIdEvent
  Bio::EnsEMBL::StableIdHistoryTree

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor;

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
our @ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

use Bio::EnsEMBL::ArchiveStableId;
use Bio::EnsEMBL::StableIdEvent;
use Bio::EnsEMBL::StableIdHistoryTree;
use Bio::EnsEMBL::Utils::Exception qw(deprecate warning throw);

use constant MAX_ROWS => 30;
use constant NUM_HIGH_SCORERS => 20;


=head2 fetch_by_stable_id

  Arg [1]     : string $stable_id
  Arg [2]     : (optional) string $type
  Example     : none
  Description : Retrives an ArchiveStableId that is the latest incarnation of
                given stable_id. If the lookup fails, attempts to check for a
                version id delimited by a period (.) and lookup again using the
                version id.
  Returntype  : Bio::EnsEMBL::ArchiveStableId or undef if not in database
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub fetch_by_stable_id {
  my $self = shift;
  my $stable_id = shift;

  my $arch_id = $self->_fetch_by_stable_id($stable_id, @_);

  # If we didn't get anything back, desperately try to see if there's
  # a version number in the stable_id
  if(!defined($arch_id) && (my $vindex = rindex($stable_id, '.'))) {
      $arch_id = $self->fetch_by_stable_id_version(substr($stable_id,0,$vindex),
						   substr($stable_id,$vindex+1),
						   @_);
  }

  return $arch_id;
}

=head2 _fetch_by_stable_id

  Arg [1]     : string $stable_id
  Arg [2]     : (optional) string $type
  Example     : none
  Description : Retrives an ArchiveStableId that is the latest incarnation of
                given stable_id. Helper function to fetch_by_stable_id, should
                not be directly called.
  Returntype  : Bio::EnsEMBL::ArchiveStableId or undef if not in database
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub _fetch_by_stable_id {
  my $self = shift;
  my $stable_id = shift;
  
  my $arch_id = Bio::EnsEMBL::ArchiveStableId->new( 
     -stable_id => $stable_id,
     -adaptor => $self
  );

  @_ ? $arch_id->type(shift) : $self->_resolve_type($arch_id);

  if ($self->lookup_current($arch_id)) {

    # stable ID is in current release
    $arch_id->version($arch_id->current_version);
    $arch_id->db_name($self->dbc->dbname);
    $arch_id->release($self->get_current_release);
    $arch_id->assembly($self->get_current_assembly);
  
  } else {

    # look for latest version of this stable id
    my $extra_sql = defined($arch_id->{'type'}) ?
      " AND sie.type = '@{[lc($arch_id->{'type'})]}'" : '';

    my $r = $self->_fetch_archive_id($stable_id, $extra_sql, $extra_sql);

    if ($r->{'new_stable_id'} and $r->{'new_stable_id'} eq $stable_id) {
      # latest event is a self event, use new_* data
      $arch_id->version($r->{'new_version'});
      $arch_id->release($r->{'new_release'});
      $arch_id->assembly($r->{'new_assembly'});
      $arch_id->db_name($r->{'new_db_name'});
    } else {
      # latest event is a deletion event (or mapping to other ID; this clause
      # is only used to cope with buggy data where deletion events are
      # missing), use old_* data
      $arch_id->version($r->{'old_version'});
      $arch_id->release($r->{'old_release'});
      $arch_id->assembly($r->{'old_assembly'});
      $arch_id->db_name($r->{'old_db_name'});
    }

    $arch_id->type(ucfirst(lc($r->{'type'})));
  }
  
  if (! defined $arch_id->db_name) {
    # couldn't find stable ID in archive or current db
    return undef;
  }

  $arch_id->is_latest(1);

  return $arch_id;
}


=head2 fetch_by_stable_id_version

  Arg [1]     : string $stable_id
  Arg [2]     : int $version
  Example     : none
  Description : Retrieve an ArchiveStableId with given version and stable ID.
  Returntype  : Bio::EnsEMBL::ArchiveStableId 
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub fetch_by_stable_id_version {
  my $self = shift;
  my $stable_id = shift;
  my $version = shift;

  my $arch_id = Bio::EnsEMBL::ArchiveStableId->new( 
     -stable_id => $stable_id,
     -version => $version,
     -adaptor => $self
  );
  
  @_ ? $arch_id->type(shift) : $self->_resolve_type($arch_id);

  if ($self->lookup_current($arch_id) && $arch_id->is_current) {

    # this version is the current one
    $arch_id->db_name($self->dbc->dbname);
    $arch_id->release($self->get_current_release);
    $arch_id->assembly($self->get_current_assembly);
  
  } else {

    # find latest release this stable ID version is found in archive
    my $extra_sql1 = qq(AND sie.old_version = "$version");
    my $extra_sql2 = qq(AND sie.new_version = "$version");
    my $r = $self->_fetch_archive_id($stable_id, $extra_sql1, $extra_sql2);

    if ($r->{'new_stable_id'} and $r->{'new_stable_id'} eq $stable_id
        and $r->{'new_version'} == $version) {
      # latest event is a self event, use new_* data
      $arch_id->release($r->{'new_release'});
      $arch_id->assembly($r->{'new_assembly'});
      $arch_id->db_name($r->{'new_db_name'});
    } else {
      # latest event is a deletion event (or mapping to other ID; this clause
      # is only used to cope with buggy data where deletion events are
      # missing), use old_* data
      $arch_id->release($r->{'old_release'});
      $arch_id->assembly($r->{'old_assembly'});
      $arch_id->db_name($r->{'old_db_name'});
    }

    $arch_id->type(ucfirst(lc($r->{'type'})));
  }
  
  if (! defined $arch_id->db_name) {
    # couldn't find stable ID version in archive or current release
    return undef;
  }

  return $arch_id;
}


=head2 fetch_by_stable_id_dbname

  Arg [1]     : string $stable_id
  Arg [2]     : string $db_name
  Example     : none
  Description : Create an ArchiveStableId from given arguments.
  Returntype  : Bio::EnsEMBL::ArchiveStableId or undef if not in database
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub fetch_by_stable_id_dbname {
  my $self = shift;
  my $stable_id = shift;
  my $db_name = shift;
  
  my $arch_id = Bio::EnsEMBL::ArchiveStableId->new( 
     -stable_id => $stable_id,
     -db_name => $db_name,
     -adaptor => $self
  );
  
  @_ ? $arch_id->type(shift) : $self->_resolve_type($arch_id);

  if ($self->lookup_current($arch_id) and $db_name eq $self->dbc->dbname) {

    # this version is the current one
    $arch_id->version($arch_id->current_version);
    $arch_id->release($self->get_current_release);
    $arch_id->assembly($self->get_current_assembly);
  
  } else {

    # find version for this dbname in the stable ID archive
    my $extra_sql = defined($arch_id->{'type'}) ?
      " AND sie.type = '@{[lc($arch_id->{'type'})]}'" : '';
    my $extra_sql1 = $extra_sql . qq( AND ms.old_db_name = "$db_name");
    my $extra_sql2 = $extra_sql . qq( AND ms.new_db_name = "$db_name");
    my $r = $self->_fetch_archive_id($stable_id, $extra_sql1, $extra_sql2);

    if ($r->{'new_stable_id'} and $r->{'new_stable_id'} eq $stable_id
        and $r->{'new_db_name'} eq $db_name) {

      # latest event is a self event, use new_* data
      $arch_id->release($r->{'new_release'});
      $arch_id->assembly($r->{'new_assembly'});
      $arch_id->version($r->{'new_version'});
    } else {
      # latest event is a deletion event (or mapping to other ID; this clause
      # is only used to cope with buggy data where deletion events are
      # missing), use old_* data
      $arch_id->release($r->{'old_release'});
      $arch_id->assembly($r->{'old_assembly'});
      $arch_id->version($r->{'old_version'});
    }

    $arch_id->type(ucfirst(lc($r->{'type'})));
  }
  
  if (! defined $arch_id->version ) {
    # couldn't find stable ID version in archive or current release
    return undef;
  }

  return $arch_id;
}

#
# Helper method to do fetch ArchiveStableId from db.
# Used by fetch_by_stable_id(), fetch_by_stable_id_version() and
# fetch_by_stable_id_dbname().
# Returns hashref as returned by DBI::sth::fetchrow_hashref
#
sub _fetch_archive_id {
  my $self = shift;
  my $stable_id = shift;
  my $extra_sql1 = shift;
  my $extra_sql2 = shift;

  # using a UNION is much faster in this query than somthing like
  # "... AND (sie.old_stable_id = ? OR sie.new_stable_id = ?)"
  my $sql = qq(
    SELECT * FROM stable_id_event sie, mapping_session ms
    WHERE sie.mapping_session_id = ms.mapping_session_id
    AND sie.old_stable_id = ?
    $extra_sql1
    UNION
    SELECT * FROM stable_id_event sie, mapping_session ms
    WHERE sie.mapping_session_id = ms.mapping_session_id
    AND sie.new_stable_id = ?
    $extra_sql2
    ORDER BY created DESC, score DESC
    LIMIT 1
  );

  my $sth = $self->prepare($sql);
  $sth->execute($stable_id,$stable_id);
  my $r = $sth->fetchrow_hashref;
  $sth->finish;

  return $r;
}  


=head2 fetch_all_by_archive_id

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId $archive_id
  Arg [2]     : String $return_type - type of ArchiveStableId to fetch
  Example     : my $arch_id = $arch_adaptor->fetch_by_stable_id('ENSG0001');
                my @archived_transcripts =
                 $arch_adaptor->fetch_all_by_archive_id($arch_id, 'Transcript');
  Description : Given a ArchiveStableId it retrieves associated ArchiveStableIds
                of specified type (e.g. retrieve transcripts for genes or vice
                versa).

                See also fetch_associated_archived() for a different approach to
                retrieve this data.
  Returntype  : listref Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : Bio::EnsEMBL::ArchiveStableId->get_all_gene_archive_ids,
                get_all_transcript_archive_ids, get_all_translation_archive_ids
  Status      : At Risk
              : under development

=cut

sub fetch_all_by_archive_id {
  my $self = shift;
  my $archive_id = shift;
  my $return_type = shift;

  my @result = ();
  my $lc_self_type = lc($archive_id->type);
  my $lc_return_type = lc($return_type);

  my $sql = qq(
    SELECT
          ga.${lc_return_type}_stable_id,
          ga.${lc_return_type}_version,
          m.old_db_name,
          m.old_release,
          m.old_assembly
    FROM  gene_archive ga, mapping_session m
    WHERE ga.${lc_self_type}_stable_id = ?
    AND   ga.${lc_self_type}_version = ?
    AND   ga.mapping_session_id = m.mapping_session_id
  );
  
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $archive_id->stable_id, SQL_VARCHAR);
  $sth->bind_param(2, $archive_id->version, SQL_SMALLINT);
  $sth->execute;
  
  my ($stable_id, $version, $db_name, $release, $assembly);
  $sth->bind_columns(\$stable_id, \$version, \$db_name, \$release, \$assembly);

  while ($sth->fetch) {
    my $new_arch_id = Bio::EnsEMBL::ArchiveStableId->new(
       -stable_id => $stable_id,
       -version => $version,
       -db_name => $db_name,
       -release => $release,
       -assembly => $assembly,
       -type => $return_type,
       -adaptor => $self
    );

    push( @result, $new_arch_id );
  }

  $sth->finish();
  return \@result;
}


=head2 fetch_associated_archived 

  Arg[1]      : Bio::EnsEMBL::ArchiveStableId $arch_id -
                the ArchiveStableId to fetch associated archived IDs for
  Example     : my ($arch_gene, $arch_tr, $arch_tl, $pep_seq) =
                  @{ $archive_adaptor->fetch_associated_archived($arch_id) };
  Description : Fetches associated archived stable IDs from the db for a given
                ArchiveStableId (version is taken into account).
  Return type : Listref of
                  ArchiveStableId archived gene
                  ArchiveStableId archived transcript
                  (optional) ArchiveStableId archived translation
                  (optional) peptide sequence
  Exceptions  : thrown on missing or wrong argument
                thrown if ArchiveStableID has no type
  Caller      : Bio::EnsEMBL::ArchiveStableId->get_all_associated_archived()
  Status      : At Risk
              : under development

=cut

sub fetch_associated_archived {
  my $self = shift;
  my $arch_id = shift;

  throw("Need a Bio::EnsEMBL::ArchiveStableId") unless ($arch_id
    and ref($arch_id) and $arch_id->isa('Bio::EnsEMBL::ArchiveStableId'));

  my $type = $arch_id->type();

  if ( !defined($type) ) {
    throw("Can't deduce ArchiveStableId type.");
  }

  $type = lc($type);

  my $sql = qq(
    SELECT  ga.gene_stable_id,
            ga.gene_version,
            ga.transcript_stable_id,
            ga.transcript_version,
            ga.translation_stable_id,
            ga.translation_version,
            pa.peptide_seq,
            ms.old_release,
            ms.old_assembly,
            ms.old_db_name
    FROM (mapping_session ms, gene_archive ga)
    LEFT JOIN peptide_archive pa
      ON ga.peptide_archive_id = pa.peptide_archive_id
    WHERE ga.mapping_session_id = ms.mapping_session_id
    AND ga.${type}_stable_id = ?
    AND ga.${type}_version = ?
  );

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $arch_id->stable_id, SQL_VARCHAR);
  $sth->bind_param(2, $arch_id->version, SQL_SMALLINT);
  $sth->execute;

  my @result = ();

  while (my $r = $sth->fetchrow_hashref) {

    my @row = ();

    # create ArchiveStableIds genes, transcripts and translations
    push @row, Bio::EnsEMBL::ArchiveStableId->new(
       -stable_id => $r->{'gene_stable_id'},
       -version => $r->{'gene_version'},
       -db_name => $r->{'old_db_name'},
       -release => $r->{'old_release'},
       -assembly => $r->{'old_assembly'},
       -type => 'Gene',
       -adaptor => $self
    );
    
    push @row, Bio::EnsEMBL::ArchiveStableId->new(
       -stable_id => $r->{'transcript_stable_id'},
       -version => $r->{'transcript_version'},
       -db_name => $r->{'old_db_name'},
       -release => $r->{'old_release'},
       -assembly => $r->{'old_assembly'},
       -type => 'Transcript',
       -adaptor => $self
    );

    if ($r->{'translation_stable_id'}) {
      push @row, Bio::EnsEMBL::ArchiveStableId->new(
         -stable_id => $r->{'translation_stable_id'},
         -version => $r->{'translation_version'},
         -db_name => $r->{'old_db_name'},
         -release => $r->{'old_release'},
         -assembly => $r->{'old_assembly'},
         -type => 'Translation',
         -adaptor => $self
      );

      # push peptide sequence onto result list
      push @row, $r->{'peptide_seq'};
    }
    
    push @result, \@row;
  }

  return \@result;
}


=head2 fetch_predecessors_by_archive_id

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId
  Example     : none
  Description : Retrieve a list of ArchiveStableIds that were mapped to the 
                given one. This method goes back only one level, to retrieve
                a full predecessor history use fetch_predecessor_history, or 
                ideally fetch_history_tree_by_stable_id for the complete
                history network.
  Returntype  : listref Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : Bio::EnsEMBL::ArchiveStableId->get_all_predecessors
  Status      : At Risk
              : under development

=cut

sub fetch_predecessors_by_archive_id {
  my $self = shift;
  my $arch_id = shift;
  
  my @result;
  
  if( ! ( defined $arch_id->stable_id() &&
	  defined $arch_id->db_name() )) {
    throw( "Need db_name for predecessor retrieval" );
  }

  my $sql = qq(
    SELECT
          sie.old_stable_id,
          sie.old_version,
          sie.type,
          m.old_db_name,
          m.old_release,
          m.old_assembly
    FROM  mapping_session m, stable_id_event sie
    WHERE sie.mapping_session_id = m.mapping_session_id
    AND   sie.new_stable_id = ?
    AND   m.new_db_name = ?	
  );

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $arch_id->stable_id, SQL_VARCHAR);
  $sth->bind_param(2, $arch_id->db_name, SQL_VARCHAR);
  $sth->execute();
  
  my ($old_stable_id, $old_version, $type, $old_db_name, $old_release, $old_assembly);
  $sth->bind_columns(\$old_stable_id, \$old_version, \$type, \$old_db_name, \$old_release, \$old_assembly);
  
  while ($sth->fetch) {
    if (defined $old_stable_id) {
      my $old_arch_id = Bio::EnsEMBL::ArchiveStableId->new( 
	 -stable_id => $old_stable_id,
	 -version => $old_version,
	 -db_name => $old_db_name,
         -release => $old_release,
         -assembly => $old_assembly,
         -type => $type,
	 -adaptor => $self
      );
      push( @result, $old_arch_id );
    }
  }
  $sth->finish();

  # if you didn't find any predecessors, there might be a gap in the
  # mapping_session history (i.e. databases in mapping_session don't chain). To
  # bridge the gap, look in the previous mapping_session for identical
  # stable_id.version
  unless (@result) {

    $sql = qq(
      SELECT
            sie.new_stable_id,
            sie.new_version,
            sie.type,
            m.new_db_name,
            m.new_release,
            m.new_assembly
      FROM  mapping_session m, stable_id_event sie
      WHERE sie.mapping_session_id = m.mapping_session_id
      AND   sie.new_stable_id = ?
      AND   m.new_db_name = ?	
    );

    $sth = $self->prepare($sql);

    my $curr_dbname = $arch_id->db_name;
    
    PREV:
    while (my $prev_dbname = $self->previous_dbname($curr_dbname)) {
    
      $sth->bind_param(1,$arch_id->stable_id, SQL_VARCHAR);
      $sth->bind_param(2,$prev_dbname, SQL_VARCHAR);
      $sth->execute();
      
      $sth->bind_columns(\$old_stable_id, \$old_version, \$type, \$old_db_name, \$old_release, \$old_assembly);
      
      while( $sth->fetch() ) {
        if (defined $old_stable_id) {
          my $old_arch_id = Bio::EnsEMBL::ArchiveStableId->new( 
             -stable_id => $old_stable_id,
             -version => $old_version,
             -db_name => $old_db_name,
             -release => $old_release,
             -assembly => $old_assembly,
             -type => $type,
             -adaptor => $self
          );
          push( @result, $old_arch_id );

          last PREV;
        }
      }

      $curr_dbname = $prev_dbname;

    }
      
    $sth->finish();
  }

  return \@result;
}


=head2 fetch_successors_by_archive_id

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId
  Example     : none
  Description : Retrieve a list of ArchiveStableIds that the given one was 
                mapped to. This method goes forward only one level, to retrieve
                a full successor history use fetch_successor_history, or 
                ideally fetch_history_tree_by_stable_id for the complete
                history network.
  Returntype  : listref Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : Bio::EnsEMBL::ArchiveStableId->get_all_successors
  Status      : At Risk
              : under development

=cut

sub fetch_successors_by_archive_id {
  my $self = shift;
  my $arch_id = shift;
  my @result;

  
  if( ! ( defined $arch_id->stable_id() &&
	  defined $arch_id->db_name() )) {
    throw( "Need db_name for successor retrieval" );
  }

  my $sql = qq(
    SELECT
          sie.new_stable_id,
          sie.new_version,
          sie.type,
          m.new_db_name,
          m.new_release,
          m.new_assembly
    FROM  mapping_session m, stable_id_event sie
    WHERE sie.mapping_session_id = m.mapping_session_id
    AND   sie.old_stable_id = ?
    AND   m.old_db_name = ?	
  );

  my $sth = $self->prepare( $sql );
  $sth->bind_param(1,$arch_id->stable_id,SQL_VARCHAR);
  $sth->bind_param(2,$arch_id->db_name,SQL_VARCHAR);
  $sth->execute();
  
  my ($new_stable_id, $new_version, $type, $new_db_name, $new_release, $new_assembly);
  $sth->bind_columns(\$new_stable_id, \$new_version, \$type, \$new_db_name, \$new_release, \$new_assembly);
  
  while( $sth->fetch() ) {
    if( defined $new_stable_id ) {
      my $new_arch_id = Bio::EnsEMBL::ArchiveStableId->new( 
	 -stable_id => $new_stable_id,
	 -version => $new_version,
	 -db_name => $new_db_name,
         -release => $new_release,
         -assembly => $new_assembly,
         -type => $type,
	 -adaptor => $self
      );
        
      push( @result, $new_arch_id );
    }
  }
  $sth->finish();
  
  # if you didn't find any successors, there might be a gap in the
  # mapping_session history (i.e. databases in mapping_session don't chain). To
  # bridge the gap, look in the next mapping_session for identical
  # stable_id.version
  unless (@result) {

    $sql = qq(
      SELECT
            sie.old_stable_id,
            sie.old_version,
            sie.type,
            m.old_db_name,
            m.old_release,
            m.old_assembly
      FROM  mapping_session m, stable_id_event sie
      WHERE sie.mapping_session_id = m.mapping_session_id
      AND   sie.old_stable_id = ?
      AND   m.old_db_name = ?	
    );

    $sth = $self->prepare($sql);

    my $curr_dbname = $arch_id->db_name;
    
    NEXTDB:
    while (my $next_dbname = $self->next_dbname($curr_dbname)) {

      $sth->bind_param(1, $arch_id->stable_id, SQL_VARCHAR);
      $sth->bind_param(2, $next_dbname, SQL_VARCHAR);
      $sth->execute();
      
      $sth->bind_columns(\$new_stable_id, \$new_version, \$type, \$new_db_name, \$new_release, \$new_assembly);
      
      while( $sth->fetch() ) {
        if (defined $new_stable_id) {
          my $new_arch_id = Bio::EnsEMBL::ArchiveStableId->new( 
             -stable_id => $new_stable_id,
             -version => $new_version,
             -db_name => $new_db_name,
             -release => $new_release,
             -assembly => $new_assembly,
             -type => $type,
             -adaptor => $self
          );
            
          push( @result, $new_arch_id );

          last NEXTDB;
        }
      }
      
      $curr_dbname = $next_dbname;

    }

    $sth->finish();
  }

  return \@result;
}



=head2 fetch_history_tree_by_stable_id

  Arg [1]     : String $stable_id - the stable ID to fetch the history tree for
  Arg [2]     : (optional) Int $num_high_scorers
                number of mappings per stable ID allowed when filtering
  Arg [3]     : (optional) Int $max_rows
                maximum number of stable IDs in history tree (used for
                filtering)
  Arg [4]     : (optional) Float $time_limit
                Optimise tree normally runs until it hits a minimised state
                but this can take a very long time. Therefore you can
                opt to bail out of the optimisation early. Specify the
                time in seconds. Floating point values are supported should you
                require sub-second limits
  Example     : my $history = $archive_adaptor->fetch_history_tree_by_stable_id(
                  'ENSG00023747897');
  Description : Returns the history tree for a given stable ID. This will
                include a network of all stable IDs it is related to. The
                method will try to return a minimal (sparse) set of nodes
                (ArchiveStableIds) and links (StableIdEvents) by removing any
                redundant entries and consolidating mapping events so that only
                changes are recorded.
  Return type : Bio::EnsEMBL::StableIdHistoryTree
  Exceptions  : thrown on missing argument
  Caller      : Bio::EnsEMBL::ArchiveStableId::get_history_tree, general
  Status      : At Risk
              : under development

=cut

sub fetch_history_tree_by_stable_id {
  my ($self, $stable_id, $num_high_scorers, $max_rows, $time_limit) = @_;

  throw("Expecting a stable ID argument.") unless $stable_id;

  $num_high_scorers ||= NUM_HIGH_SCORERS;
  $max_rows ||= MAX_ROWS;

  # using a UNION is much faster in this query than somthing like
  # "... AND (sie.old_stable_id = ?) OR (sie.new_stable_id = ?)"
  #
  # SQLite uses the fully qualified column name as the key in
  # fetchrow_hashref() when there's a UNION, hence the need to
  # avoid table names qualifiers in the column lists.
  my $sql = qq(
    SELECT old_stable_id, old_version,
           old_db_name, old_release, old_assembly,
           new_stable_id, new_version,
           new_db_name, new_release, new_assembly,
           type, score
    FROM stable_id_event sie, mapping_session ms
    WHERE sie.mapping_session_id = ms.mapping_session_id
    AND sie.old_stable_id = ?
    UNION
    SELECT old_stable_id, old_version,
           old_db_name, old_release, old_assembly,
           new_stable_id, new_version,
           new_db_name, new_release, new_assembly,
           type, score
    FROM stable_id_event sie, mapping_session ms
    WHERE sie.mapping_session_id = ms.mapping_session_id
    AND sie.new_stable_id = ?
  );
  
  my $sth = $self->prepare($sql);

  my $history = Bio::EnsEMBL::StableIdHistoryTree->new(
      -CURRENT_DBNAME => $self->dbc->dbname,
      -CURRENT_RELEASE => $self->get_current_release,
      -CURRENT_ASSEMBLY => $self->get_current_assembly,
  );

  # remember stable IDs you need to do and those that are done. Initialise the
  # former hash with the focus stable ID
  my %do = ($stable_id => 1);
  my %done;

  # while we got someting to do
  while (my ($id) = keys(%do)) {

    # if we already have more than MAX_ROWS stable IDs in this tree, we can't
    # build the full tree. Return undef.
    if (scalar(keys(%done)) > $max_rows) {
      # warning("Too many related stable IDs (".scalar(keys(%done)).") to draw a history tree.");
      $history->is_incomplete(1);
      $sth->finish;
      last;
    }

    # mark this stable ID as done
    delete $do{$id};
    $done{$id} = 1;

    # fetch all stable IDs related to this one from the database
    $sth->bind_param(1, $id, SQL_VARCHAR);
    $sth->bind_param(2, $id, SQL_VARCHAR);
    $sth->execute;

    my @events;

    while (my $r = $sth->fetchrow_hashref) {
      
      #
      # create old and new ArchiveStableIds and a StableIdEvent to link them
      # add all of these to the history tree
      #
      my ($old_id, $new_id);

      if ($r->{'old_stable_id'}) {
        $old_id = Bio::EnsEMBL::ArchiveStableId->new(
          -stable_id => $r->{'old_stable_id'},
          -version => $r->{'old_version'},
          -db_name => $r->{'old_db_name'},
          -release => $r->{'old_release'},
          -assembly => $r->{'old_assembly'},
          -type => $r->{'type'},
          -adaptor => $self
        );
      }
       
      if ($r->{'new_stable_id'}) {
        $new_id = Bio::EnsEMBL::ArchiveStableId->new(
          -stable_id => $r->{'new_stable_id'},
          -version => $r->{'new_version'},
          -db_name => $r->{'new_db_name'},
          -release => $r->{'new_release'},
          -assembly => $r->{'new_assembly'},
          -type => $r->{'type'},
          -adaptor => $self
        );
      }

      my $event = Bio::EnsEMBL::StableIdEvent->new(
        -old_id => $old_id,
        -new_id => $new_id,
        -score => $r->{'score'}
      );

      push @events, $event;

    }

    # filter out low-scoring events; the number of highest scoring events
    # returned is defined by NUM_HIGH_SCORERS
    my @others;

    foreach my $event (@events) {
      
      my $old_id = $event->old_ArchiveStableId;
      my $new_id = $event->new_ArchiveStableId;
      
      # creation, deletion and mapping-to-self events are added to the history
      # tree directly
      if (!$old_id || !$new_id || ($old_id->stable_id eq $new_id->stable_id)) {
        $history->add_StableIdEvents($event);
      } else {
        push @others, $event;
      }
      
    }

    #if (scalar(@others) > $num_high_scorers) {
    #  warn "Filtering ".(scalar(@others) - $num_high_scorers).
    #    " low-scoring events.\n";
    #}

    my $k = 0;
    foreach my $event (sort { $b->score <=> $a->score } @others) {
      $history->add_StableIdEvents($event);
      
      # mark stable IDs as todo if appropriate
      $do{$event->old_ArchiveStableId->stable_id} = 1
        unless $done{$event->old_ArchiveStableId->stable_id};
      $do{$event->new_ArchiveStableId->stable_id} = 1
        unless $done{$event->new_ArchiveStableId->stable_id};
      
      last if (++$k == $num_high_scorers);
    }
    
  }

  $sth->finish;
  
  # try to consolidate the tree (remove redundant nodes, bridge gaps)
  $history->consolidate_tree;

  # now add ArchiveStableIds for current Ids not found in the archive
  $self->add_all_current_to_history($history);
  
  # calculate grid coordinates for the sorted tree; this will also try to
  # untangle the tree
  $history->calculate_coords($time_limit);
  
  return $history;
}


=head2 add_all_current_to_history 

  Arg[1]      : Bio::EnsEMBL::StableIdHistoryTree $history -
                the StableIdHistoryTree object to add the current IDs to
  Description : This method adds the current versions of all stable IDs found
                in a StableIdHistoryTree object to the tree, by creating
                appropriate Events for the stable IDs found in the *_stable_id
                tables. This is a helper method for
                fetch_history_tree_by_stable_id(), see there for more
                documentation.
  Return type : none (passed-in object is manipulated)
  Exceptions  : thrown on missing or wrong argument
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub add_all_current_to_history {
  my $self = shift;
  my $history = shift;

  unless ($history and $history->isa('Bio::EnsEMBL::StableIdHistoryTree')) {
    throw("Need a Bio::EnsEMBL::StableIdHistoryTree.");
  }

  my @ids = @{ $history->get_unique_stable_ids };
  my $id_string = join("', '", @ids);
  
  my $tmp_id = Bio::EnsEMBL::ArchiveStableId->new(-stable_id => $ids[0]);
  my $type = lc($self->_resolve_type($tmp_id));
  return unless ($type);

  # get current stable IDs from db
  my $sql = qq(
    SELECT stable_id, version FROM ${type}
    WHERE stable_id IN ('$id_string')
  );
  my $sth = $self->prepare($sql);
  $sth->execute;

  while (my ($stable_id, $version) = $sth->fetchrow_array) {

    my $new_id = Bio::EnsEMBL::ArchiveStableId->new(
      -stable_id => $stable_id,
      -version => $version,
      -current_version => $version,
      -db_name => $self->dbc->dbname,
      -release => $self->get_current_release,
      -assembly => $self->get_current_assembly,
      -type => $type,
      -adaptor => $self
    );

    my $event = $history->get_latest_StableIdEvent($new_id);
    next unless ($event);

    if ($event->old_ArchiveStableId and
        $event->old_ArchiveStableId->stable_id eq $stable_id) {
      
      # latest event was a self event
      # update it with current stable ID and add to tree
      $event->new_ArchiveStableId($new_id);

    } else {

      # latest event was a non-self event
      # create a new event where the old_id is the new_id from latest
      my $new_event = Bio::EnsEMBL::StableIdEvent->new(
        -old_id => $event->new_ArchiveStableId,
        -new_id => $new_id,
        -score => $event->score,
      );
      $history->add_StableIdEvents($new_event);
    }
    
  }

  # refresh node cache
  $history->flush_ArchiveStableIds;
  $history->add_ArchiveStableIds_for_events;
}


=head2 fetch_successor_history

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId $arch_id
  Example     : none
  Description : Gives back a list of archive stable ids which are successors in
                the stable_id_event tree of the given stable_id. Might well be
                empty.
                
                This method isn't deprecated, but in most cases you will rather
                want to use fetch_history_tree_by_stable_id().
  Returntype  : listref Bio::EnsEMBL::ArchiveStableId
                Since every ArchiveStableId knows about it's successors, this is
                a linked tree.
  Exceptions  : none
  Caller      : webcode for archive
  Status      : At Risk
              : under development

=cut

sub fetch_successor_history {
  my $self = shift;
  my $arch_id = shift;

  my $current_db_name = $self->list_dbnames->[0];
  my $dbname = $arch_id->db_name;

  if ($dbname eq $current_db_name) {
    return [$arch_id];
  }

  my $old = [];
  my @result = ();
  
  push @$old, $arch_id;

  while ($dbname ne $current_db_name) {
    my $new = [];
    while (my $asi = (shift @$old)) {
      push @$new, @{ $asi->get_all_successors };
    }

    if (@$new) {
      $dbname = $new->[0]->db_name;
    } else {
      last;
    }

    # filter duplicates
    my %unique = map { join(":", $_->stable_id, $_->version, $_->release) =>
      $_ } @$new;
    @$new = values %unique;
    
    @$old = @$new;
    push @result, @$new;
  }

  return \@result;
}


=head2 fetch_predecessor_history

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId $arch_id
  Example     : none
  Description : Gives back a list of archive stable ids which are predecessors
                in the stable_id_event tree of the given stable_id. Might well
                be empty.
                
                This method isn't deprecated, but in most cases you will rather
                want to use fetch_history_tree_by_stable_id().
  Returntype  : listref Bio::EnsEMBL::ArchiveStableId
                Since every ArchiveStableId knows about it's successors, this is
                a linked tree.
  Exceptions  : none
  Caller      : webcode for archive
  Status      : At Risk
              : under development

=cut

sub fetch_predecessor_history {
  my $self = shift;
  my $arch_id = shift;

  my $oldest_db_name = $self->list_dbnames->[-1];
  my $dbname = $arch_id->db_name;

  if ($dbname eq $oldest_db_name) {
    return [$arch_id];
  }

  my $old = [];
  my @result = ();

  push @$old, $arch_id;

  while ($dbname ne $oldest_db_name) {
    my $new = [];
    while (my $asi = (shift @$old)) {
      push @$new, @{ $asi->get_all_predecessors };
    }

    if( @$new ) {
      $dbname = $new->[0]->db_name;
    } else {
      last;
    }
    
    # filter duplicates
    my %unique = map { join(":", $_->stable_id, $_->version, $_->release) =>
      $_ } @$new;
    @$new = values %unique;
    
    @$old = @$new;
    push @result, @$new;
  }

  return \@result;
}


=head2 fetch_stable_id_event

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId $arch_id
  Arg [2]     : stable_id
  Example     : my $archive = $archive_stable_id_adaptor->fetch_by_stable_id($id);
                my $event = $archive_stable_id_adaptor($archive, $id2);
  Description : Gives back the event that links an archive stable id
                to a specific stable id
                
  Returntype  : Bio::EnsEMBL::StableIdEvent
                Undef if no event was found
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut


sub fetch_stable_id_event {
  my $self = shift;
  my $arch_id = shift;
  my $stable_id = shift;

  my $event;

  my $sql = qq(
    SELECT sie.old_stable_id, sie.old_version, sie.new_stable_id, sie.new_version, sie.type, sie.score,
           ms.old_db_name, ms.new_db_name, ms.old_release, ms.new_release, ms.old_assembly, ms.new_assembly
    FROM stable_id_event sie, mapping_session ms
    WHERE ms.mapping_session_id = sie.mapping_session_id
    AND (old_stable_id = ? AND ms.old_db_name = ? AND old_release = ? AND old_assembly = ? AND new_stable_id = ?)
    OR (new_stable_id = ? AND ms.new_db_name = ? AND new_release = ? AND new_assembly = ? AND old_stable_id = ?)
  );

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $arch_id->stable_id, SQL_VARCHAR);
  $sth->bind_param(2, $arch_id->db_name, SQL_VARCHAR);
  $sth->bind_param(3, $arch_id->release, SQL_INTEGER);
  $sth->bind_param(4, $arch_id->assembly, SQL_VARCHAR);
  $sth->bind_param(5, $stable_id, SQL_VARCHAR);
  $sth->bind_param(6, $arch_id->stable_id, SQL_VARCHAR);
  $sth->bind_param(7, $arch_id->db_name, SQL_VARCHAR);
  $sth->bind_param(8, $arch_id->release, SQL_INTEGER);
  $sth->bind_param(9, $arch_id->assembly, SQL_VARCHAR);
  $sth->bind_param(10, $stable_id, SQL_VARCHAR);
  $sth->execute();

  my ($old_stable_id, $old_version, $new_stable_id, $new_version, $type, $score);
  my ($old_db_name, $new_db_name, $old_release, $new_release, $old_assembly, $new_assembly);
  $sth->bind_columns(\$old_stable_id, \$old_version, \$new_stable_id, \$new_version, \$type, \$score,
                     \$old_db_name, \$new_db_name, \$old_release, \$new_release, \$old_assembly, \$new_assembly);

  while ($sth->fetch) {
    if ($new_stable_id eq $stable_id) {

      my $alt_id = Bio::EnsEMBL::ArchiveStableId->new(
          -stable_id => $new_stable_id,
          -version => $new_version,
          -db_name => $new_db_name,
          -release => $new_release,
          -assembly => $new_assembly,
          -type => $type,
          -adaptor => $self
      );

      $event = Bio::EnsEMBL::StableIdEvent->new(
        -old_id => $arch_id,
        -new_id => $alt_id,
        -score => $score
      );

    } elsif ($old_stable_id eq $stable_id) {

      my $alt_id = Bio::EnsEMBL::ArchiveStableId->new(
          -stable_id => $old_stable_id,
          -version => $old_version,
          -db_name => $old_db_name,
          -release => $old_release,
          -assembly => $old_assembly,
          -type => $type,
          -adaptor => $self
      );

      $event = Bio::EnsEMBL::StableIdEvent->new(
        -old_id => $alt_id,
        -new_id => $arch_id,
        -score => $score
      );

    }
  }
  $sth->finish();

  return $event;
}


=head2 list_dbnames

  Args        : none
  Example     : none
  Description : A list of available database names from the latest (current) to
                the oldest (ordered).
  Returntype  : listref of strings
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub list_dbnames {
  my $self = shift;
  
  if( ! defined $self->{'dbnames'} ) {

    my $sql = qq(
      SELECT old_db_name, new_db_name
      FROM mapping_session
      ORDER BY created DESC
    );
    my $sth = $self->prepare( $sql );
    $sth->execute();
    my ( $old_db_name, $new_db_name );
    
    my @dbnames = ();
    my %seen;

    $sth->bind_columns( \$old_db_name, \$new_db_name );

    while( $sth->fetch() ) {
      # this code now can deal with non-chaining mapping sessions
      push(@{ $self->{'dbnames'} }, $new_db_name) unless ($seen{$new_db_name});
      $seen{$new_db_name} = 1;

      push(@{ $self->{'dbnames'} }, $old_db_name) unless ($seen{$old_db_name});
      $seen{$old_db_name} = 1;
    }

    $sth->finish();
    
  }

  return $self->{'dbnames'};
}


=head2 previous_dbname

  Arg[1]      : String $dbname - focus db name
  Example     : my $prev_db = $self->previous_dbname($curr_db);
  Description : Returns the name of the next oldest database which has mapping
                session information.
  Return type : String (or undef if not available)
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub previous_dbname {
  my $self = shift;
  my $dbname = shift;

  my $curr_idx = $self->_dbname_index($dbname);
  my @dbnames = @{ $self->list_dbnames };

  if ($curr_idx == @dbnames) {
    # this is the oldest dbname, so no previous one available
    return undef;
  } else {
    return $dbnames[$curr_idx+1];
  }
}


=head2 next_dbname

  Arg[1]      : String $dbname - focus db name
  Example     : my $prev_db = $self->next_dbname($curr_db);
  Description : Returns the name of the next newest database which has mapping
                session information.
  Return type : String (or undef if not available)
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub next_dbname {
  my $self = shift;
  my $dbname = shift;

  my $curr_idx = $self->_dbname_index($dbname);
  my @dbnames = @{ $self->list_dbnames };

  if ($curr_idx == 0) {
    # this is the latest dbname, so no next one available
    return undef;
  } else {
    return $dbnames[$curr_idx-1];
  }
}


#
# helper method to return the array index of a database in the ordered list of
# available databases (as returned by list_dbnames()
#
sub _dbname_index {
  my $self = shift;
  my $dbname = shift;

  my @dbnames = @{ $self->list_dbnames };

  for (my $i = 0; $i < @dbnames; $i++) {
    if ($dbnames[$i] eq $dbname) {
      return $i;
    }
  }
}


=head2 get_peptide

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId $arch_id
  Example     : none
  Description : Retrieves the peptide string for given ArchiveStableId. If its
                not a peptide or not in the database returns undef.
  Returntype  : string or undef
  Exceptions  : none
  Caller      : Bio::EnsEMBL::ArchiveStableId->get_peptide, general
  Status      : At Risk
              : under development

=cut

sub get_peptide {
  my $self    = shift;
  my $arch_id = shift;

  if ( lc( $arch_id->type() ) ne 'translation' ) {
    return undef;
  }

  my $sql = qq(
    SELECT pa.peptide_seq
      FROM peptide_archive pa, gene_archive ga
     WHERE ga.translation_stable_id = ?
       AND ga.translation_version = ?
       AND ga.peptide_archive_id = pa.peptide_archive_id
  );

  my $sth = $self->prepare($sql);
  $sth->bind_param( 1, $arch_id->stable_id, SQL_VARCHAR );
  $sth->bind_param( 2, $arch_id->version,   SQL_SMALLINT );
  $sth->execute();

  my ($peptide_seq) = $sth->fetchrow_array();
  $sth->finish();

  return $peptide_seq;
} ## end sub get_peptide


=head2 get_current_release 

  Example     : my $current_release = $archive_adaptor->get_current_release;
  Description : Returns the current release number (as found in the meta table).
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_current_release {
  my $self = shift;

  unless ($self->{'current_release'}) {
    my $mca = $self->db->get_MetaContainer;
    my ($release) = @{ $mca->list_value_by_key('schema_version') };
    $self->{'current_release'} = $release;
  }

  return $self->{'current_release'};
}


=head2 get_current_assembly

  Example     : my $current_assembly = $archive_adaptor->get_current_assembly;
  Description : Returns the current assembly version (as found in the meta
                table).
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_current_assembly {
  my $self = shift;

  unless ($self->{'current_assembly'}) {
    my $mca = $self->db->get_MetaContainer;
    my ($assembly) = @{ $mca->list_value_by_key('assembly.default') };
    $self->{'current_assembly'} = $assembly;
  }

  return $self->{'current_assembly'};
}


=head2 lookup_current

  Arg[1]      : Bio::EnsEMBL::ArchiveStableId $arch_id -
                the stalbe ID to find the current version for
  Example     : if ($self->lookup_version($arch_id) {
                  $arch_id->version($arch_id->current_version);
                  $arch_id->db_name($self->dbc->dbname);
  Description : Look in [gene|transcript|translation]_stable_id if you can find
                a current version for this stable ID. Set
                ArchiveStableId->current_version if found.
  Return type : Boolean (TRUE if current version found, else FALSE)
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub lookup_current {
  my $self    = shift;
  my $arch_id = shift;

  my $type = lc( $arch_id->type );

  unless ($type) {
    warning("Can't lookup current version without a type.");
    return 0;
  }

  my $sql = qq(
    SELECT version FROM ${type}
    WHERE stable_id = ?
  );
  my $sth = $self->prepare($sql);
  $sth->execute( $arch_id->stable_id );
  my ($version) = $sth->fetchrow_array;
  $sth->finish;

  if ($version) {
    $arch_id->current_version($version);
    return 1;
  }

  # didn't find a current version
  return 0;
} ## end sub lookup_current


# infer type from stable ID format
sub _resolve_type {
  my $self = shift;
  my $arch_id = shift;
  
  my $stable_id = $arch_id->stable_id();
  my $id_type;

  # first, try to infer type from stable ID format
  #
  # Anopheles IDs
  if ($stable_id =~ /^AGAP.*/) {
    if ($stable_id =~ /.*-RA/) {
      $id_type = "Transcript";
    } elsif ($stable_id =~ /.*-PA/) {
      $id_type = "Translation";
    } else {
      $id_type = "Gene";
    }

  # standard Ensembl IDs
  } elsif ($stable_id =~ /.*G\d+(\.\d+)?$/) {
    $id_type = "Gene";
  } elsif ($stable_id =~ /.*T\d+(\.\d+)?$/) { 
    $id_type = "Transcript";
  } elsif ($stable_id =~ /.*P\d+(\.\d+)?$/) { 
    $id_type = "Translation";
  } elsif ($stable_id =~ /.*E\d+(\.\d+)?$/) { 
    $id_type = "Exon";

  # if guessing fails, look in db
  } else {
    my $sql = qq(
      SELECT type from stable_id_event
      WHERE old_stable_id = ?
      OR new_stable_id = ?
    );
    my $sth = $self->prepare($sql);
    $sth->execute($stable_id, $stable_id);
    ($id_type) = $sth->fetchrow_array;
    $sth->finish;
  }

  warning("Couldn't resolve stable ID type.") unless ($id_type);
  
  $arch_id->type($id_type);
}


1;


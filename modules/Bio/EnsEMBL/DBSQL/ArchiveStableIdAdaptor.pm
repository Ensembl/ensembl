package Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor;

=head1 NAME

Bio::EnsEMBL::ArchiveStableIdAdaptor

=head1 SYNOPSIS

my $reg = "Bio::EnsEMBL::Registry";
my $archiveStableIdAdaptor =
  $reg->get_adaptor('human', 'core', 'ArchiveStableId');
  
my $arch_id = $archiveStableIdAdaptor->fetch_by_stable_id("ENSG00000068990");
my @history = @{ $archiveStableIdAdaptor->fetch_stable_id_history($arch_id) };

foreach my $a (@history) {
  print "Stable ID: ".$a->stable_id.".".$a->version."\n";
  print "Release: ".$a->release." (".$a->assembly.", ".$a->db_name.")\n");
}

=head1 DESCRIPTION

ArchiveStableIdAdaptor does all SQL to create ArchiveStableIds and works of 

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
    fetch_stable_id_history
    fetch_predecessor_history
    fetch_successor_history
    get_peptide
    list_dbnames

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Ensembl core API team
Currently maintained by Patrick Meidl <meidl@ebi.ac.uk>

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
our @ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

use Bio::EnsEMBL::ArchiveStableId;
use Bio::EnsEMBL::StableIdEvent;
use Bio::EnsEMBL::StableIdHistoryTree;
use Bio::EnsEMBL::Utils::Exception qw(deprecate warning);

use Data::Dumper;


=head2 fetch_by_stable_id

  Arg [1]     : string $stable_id
  Example     : none
  Description : retrives an ArchiveStableId that is the latest incarnation of
                given stable_id.
  Returntype  : Bio::EnsEMBL::ArchiveStableId or undef if not in database
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub fetch_by_stable_id {
  my $self = shift;
  my $stable_id = shift;
  
  my $arch_id = Bio::EnsEMBL::ArchiveStableId->new
    ( 
     -adaptor => $self,
     -stable_id => $stable_id
    );

  @_ ? $arch_id->type( shift ) : _resolve_type( $arch_id );

  if( ! $self->_lookup_version( $arch_id ) ) {
    return undef;
  }

  return $arch_id;
}


=head2 fetch_by_stable_id_version

  Arg [1]     : string $stable_id
  Arg [2]     : int $version
  Example     : none
  Description : Retrieve an archiveStableId with given version and stable ID.
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

  my $arch_id = Bio::EnsEMBL::ArchiveStableId->new
    ( 
     -adaptor => $self,
     -version => $version,
     -stable_id => $stable_id
    );
  
  @_ ? $arch_id->type( shift ) : _resolve_type( $arch_id );

  # find latest release this stable ID is found
  my $sql_tmp = qq(
    SELECT
          m.new_db_name,
          m.new_release,
          m.new_assembly,
          sie.score
    FROM  stable_id_event sie, mapping_session m
    WHERE sie.mapping_session_id = m.mapping_session_id
    AND   sie.new_stable_id = "$stable_id"
    AND   sie.new_version = $version
    AND   m.new_db_name <> 'LATEST'
    ORDER BY m.created DESC
  );
  my $sql = $self->dbc->add_limit_clause($sql_tmp, 1);

  my $sth = $self->prepare($sql);
  $sth->execute();
  my ($db_name, $release, $assembly, $score) = $sth->fetchrow_array();
  $sth->finish();
  
  # you might have missed a stable ID that was deleted in the very first
  # stable ID mapping for this species, so go back and try again
  if( ! defined $db_name ) {
    my $sql_tmp = qq(
      SELECT
            m.old_db_name,
            m.old_release,
            m.old_assembly,
            sie.score
      FROM  stable_id_event sie, mapping_session m
      WHERE sie.mapping_session_id = m.mapping_session_id
      AND   sie.old_stable_id = "$stable_id"
      AND   sie.old_version = $version
      AND   m.old_db_name <> 'ALL'
      ORDER BY m.created DESC
    );
    $sql = $self->dbc->add_limit_clause($sql_tmp, 1);

    $sth = $self->prepare($sql);
    $sth->execute();
    ($db_name, $release, $assembly, $score) = $sth->fetchrow_array();
    $sth->finish();
  }
  
  if( ! defined $db_name ) {
    # couldn't find stable ID version in archive
    return undef;
  } else {
    $arch_id->db_name($db_name);
    $arch_id->release($release);
    $arch_id->assembly($assembly);
    $arch_id->score($score);
  }

  return $arch_id;
}


=head2 fetch_by_stable_id_dbname

  Arg [1]     : string $stable_id
  Arg [2]     : string $db_name
  Example     : none
  Description : create an ArchiveStableId from given arguments.
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
  
  my $arch_id = Bio::EnsEMBL::ArchiveStableId->new
    ( 
     -adaptor => $self,
     -db_name => $db_name,
     -stable_id => $stable_id
    );
  
  @_ ? $arch_id->type( shift ) : _resolve_type( $arch_id );

  if( ! $self->_lookup_version( $arch_id ) ) {
    return undef;
  }

  return $arch_id;
}


=head2 fetch_all_by_archive_id

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId $archive_id
  Arg [2]     : String $return_type - type of ArchiveStableId to fetch
  Example     : my $arch_id = $arch_adaptor->fetch_by_stable_id('ENSG0001');
                my @archived_transcripts = $arch_adaptor->fetch_all_by_archive_id($arch_id, 'Transcript');
  Description : Given a ArchiveStableId it retrieves associated ArchiveStableIds
                of specified type (e.g. retrieve transcripts for genes or vice
                versa).
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
    my $new_arch_id = Bio::EnsEMBL::ArchiveStableId->new
      (
       -version => $version,
       -adaptor => $self,
       -stable_id => $stable_id,
       -type => $return_type,
       -db_name => $db_name,
       -release => $release,
       -assembly => $assembly
      );

    push( @result, $new_arch_id );
  }

  $sth->finish();
  return \@result;
}


=head2 fetch_predecessors_by_archive_id

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId
  Example     : none
  Description : Retrieve a list of ArchiveStableIds that were mapped to the 
                given one. 
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
    $self->throw( "Need db_name for predecessor retrieval" );
  }

  my $sql = qq(
    SELECT
          sie.old_stable_id,
          sie.old_version,
          sie.score,
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
  
  my ($old_stable_id, $old_version, $score, $old_db_name, $old_release, $old_assembly);
  $sth->bind_columns(\$old_stable_id, \$old_version, \$score, \$old_db_name, \$old_release, \$old_assembly);
  
  while( $sth->fetch() ) {
    if( defined $old_stable_id ) {
      my $old_arch_id = Bio::EnsEMBL::ArchiveStableId->new
	( 
	 -version => $old_version,
	 -stable_id => $old_stable_id,
	 -db_name => $old_db_name,
         -release => $old_release,
         -assembly => $old_assembly,
	 -adaptor => $self,
         -score => $score
	);
      _resolve_type( $old_arch_id );
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
            sie.score,
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
      
      $sth->bind_columns(\$old_stable_id, \$old_version, \$score, \$old_db_name, \$old_release, \$old_assembly);
      
      while( $sth->fetch() ) {
        if (defined $old_stable_id) {
          my $old_arch_id = Bio::EnsEMBL::ArchiveStableId->new
            ( 
             -version => $old_version,
             -stable_id => $old_stable_id,
             -db_name => $old_db_name,
             -release => $old_release,
             -assembly => $old_assembly,
             -adaptor => $self,
             -score => $score
            );
          _resolve_type( $old_arch_id );
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
                a full successor history use fetch_successor_history().
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
    $self->throw( "Need db_name for successor retrieval" );
  }

  my $sql = qq(
    SELECT
          sie.new_stable_id,
          sie.new_version,
          sie.score,
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
  
  my ($new_stable_id, $new_version, $score, $new_db_name, $new_release, $new_assembly);
  $sth->bind_columns(\$new_stable_id, \$new_version, \$score, \$new_db_name, \$new_release, \$new_assembly);
  
  while( $sth->fetch() ) {
    if( defined $new_stable_id ) {
      my $new_arch_id = Bio::EnsEMBL::ArchiveStableId->new
	( 
	 -version => $new_version,
	 -stable_id => $new_stable_id,
	 -db_name => $new_db_name,
         -release => $new_release,
         -assembly => $new_assembly,
	 -adaptor => $self,
         -score => $score
	);
        
      _resolve_type($new_arch_id);
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
            sie.score,
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
      
      $sth->bind_columns(\$new_stable_id, \$new_version, \$score, \$new_db_name, \$new_release, \$new_assembly);
      
      while( $sth->fetch() ) {
        if (defined $new_stable_id) {
          my $new_arch_id = Bio::EnsEMBL::ArchiveStableId->new
            ( 
             -version => $new_version,
             -stable_id => $new_stable_id,
             -db_name => $new_db_name,
             -release => $new_release,
             -assembly => $new_assembly,
             -adaptor => $self,
             -score => $score
            );
            
          _resolve_type($new_arch_id);
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

  Arg[1]      : String $stable_id - the stable ID to fetch the history tree for
  Example     : my $history = $archive_adaptor->fetch_history_tree_by_stable_id(
                  'ENSG00023747897');
  Description : Returns the history tree for a given stable ID. This will
                include a network of all stable IDs this ID is related to. The
                method will try to return a minimal (sparse) set of nodes
                (ArchiveStableIds) and links (StableIdEvents) by removing any
                redundant entries and consolidating mapping events so that only
                changes are recorded.
  Return type : Bio::EnsEMBL::StableIdHistoryTree
  Exceptions  : 
  Caller      : Bio::EnsEMBL::ArchiveStableId::get_history_tree, general
  Status      : At Risk
              : under development

=cut

sub fetch_history_tree_by_stable_id {
  my ($self, $stable_id) = @_;

  throw("Expecting a stable ID argument.") unless $stable_id;

  # using a UNION is much faster in this query than somthing like
  # "... AND (sie.old_stable_id = ?) OR (sie.new_stable_id = ?)"
  my $sql = qq(
    SELECT sie.old_stable_id, sie.old_version,
           ms.old_db_name, ms.old_release, ms.old_assembly,
           sie.new_stable_id, sie.new_version,
           ms.new_db_name, ms.new_release, ms.new_assembly,
           sie.type, sie.score
    FROM stable_id_event sie, mapping_session ms
    WHERE sie.mapping_session_id = ms.mapping_session_id
    AND sie.old_stable_id = ?
    UNION
    SELECT sie.old_stable_id, sie.old_version,
           ms.old_db_name, ms.old_release, ms.old_assembly,
           sie.new_stable_id, sie.new_version,
           ms.new_db_name, ms.new_release, ms.new_assembly,
           sie.type, sie.score
    FROM stable_id_event sie, mapping_session ms
    WHERE sie.mapping_session_id = ms.mapping_session_id
    AND sie.new_stable_id = ?
  );
  
  my $sth = $self->prepare($sql);

  my $history = Bio::EnsEMBL::StableIdHistoryTree->new;

  # remember stable IDs you need to do and those that are done. Initialise the
  # former hash with the focus stable ID
  my %do = ($stable_id => 1);
  my %done;

  # while we got someting to do
  while (my ($id) = keys(%do)) {

    warn "$id\n";

    # mark this stable ID as done
    delete $do{$id};
    $done{$id} = 1;

    # fetch all stable IDs related to this one from the database
    $sth->bind_param(1, $id, SQL_VARCHAR);
    $sth->bind_param(2, $id, SQL_VARCHAR);
    $sth->execute;

    while (my $r = $sth->fetchrow_hashref) {
      
      #
      # create old and new ArchiveStableIds and a StableIdEvent to link them
      # add all of these to the history tree
      #
      my ($old_id, $new_id);

      Data::Dumper::Dumper($r);
      
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
        
        # add to history tree
        $history->add_ArchiveStableIds($old_id);
        
        # mark stable IDs as todo if appropriate
        $do{$old_id->stable_id} = 1 unless $done{$old_id->stable_id};
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

        # add to history tree
        $history->add_ArchiveStableIds($new_id);
        
        # mark stable IDs as todo if appropriate
        $do{$new_id->stable_id} = 1 unless $done{$new_id->stable_id};
      }

      my $event = Bio::EnsEMBL::StableIdEvent->new(
        -old_id => $old_id,
        -new_id => $new_id,
        -score => $r->{'score'}
      );

      # add to history tree
      $history->add_StableIdEvents($event);

    }
  }

  $sth->finish;
  
  #
  # now try to consolidate the tree
  #
  # this will remove any nodes where there were no changes, connect the
  # affected links, and also create links for implicit mappings (i.e. bridge
  # gaps in the history)
  #
  # [todo]

  
  # calculate coordinates for the sorted tree
  $history->calculate_simple_coords;
  
  return $history;
}


=head2 fetch_archive_id_history

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId $arch_id
  Example     : none
  Description : Gives back a list of archive stable ids which are successors or
                predecessors in the stable_id_event tree of the given
                stable_id. Might well be empty. This is not the complete network
                this stable id belongs to, but rather branches out from this id
                only.
  Returntype  : listref of Bio::EnsEMBL::ArchiveStableId
                Since every ArchiveStableId knows about it's successors, this is
                a linked tree.
  Exceptions  : none
  Caller      : webcode for archive
  Status      : At Risk
              : under development

=cut

sub fetch_archive_id_history {
  my $self = shift;
  my $arch_id = shift;

  my @result = (
    $arch_id,
    @{ $self->fetch_predecessor_history($arch_id) },
    @{ $self->fetch_successor_history($arch_id) }
  );

  # filter duplicates
  my %unique = map { join(":", $_->stable_id, $_->version, $_->release) => $_ }
    @result;
  @result = values %unique;

  return \@result;
}


=head2 fetch_successor_history

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId $arch_id
  Example     : none
  Description : Gives back a list of archive stable ids which are successors in
                the stable_id_event tree of the given stable_id. Might well be
                empty.
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

      if( $old_db_name eq "ALL" ) {
	next;
      }

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
  my $self = shift;
  my $arch_id = shift;

  if( ! $arch_id->{'version'} ) {
    if(  ! $ self->_lookup_version( $arch_id )) {
      return undef;
    }
  }

  
  my $sql = qq(
    SELECT pa.peptide_seq
      FROM peptide_archive pa, gene_archive ga
     WHERE ga.translation_stable_id = ?
       AND ga.translation_version = ?
       AND ga.peptide_archive_id = pa.peptide_archive_id
  );


  my $sth = $self->prepare( $sql );
  $sth->bind_param(1,$arch_id->stable_id, SQL_VARCHAR);
  $sth->bind_param(2,$arch_id->version, SQL_SMALLINT);
  $sth->execute();
  
  my ( $peptide_seq ) = $sth->fetchrow_array();
  $sth->finish();

  return $peptide_seq;
}


# given an ArchiveStableId that's missing the version but has the db_name
# this should fill in the version
# return true or false 
sub _lookup_version {
  my $self = shift;
  my $arch_id = shift;

  if( defined $arch_id->{'version'} ) {
    return 1;
  }

  my $sql;
  my $EXTRA_SQL = defined($arch_id->{'type'}) ?
    " AND sie.type = '@{[lc($arch_id->{'type'})]}'" : '';

  if( ! defined $arch_id->{'db_name'} ) {
    # latest version of this stable id
    my $sql_tmp = qq(
      SELECT
            m.new_db_name,
            m.new_release,
            m.new_assembly,
            sie.new_version,
            sie.score
      FROM  stable_id_event sie, mapping_session m
      WHERE sie.mapping_session_id = m.mapping_session_id
      AND   new_stable_id = "@{[$arch_id->stable_id]}"
      AND   m.new_db_name <> 'LATEST'
      $EXTRA_SQL
      ORDER BY m.created DESC
    );
    $sql = $self->dbc->add_limit_clause($sql_tmp, 1);
  } else {
    $sql = qq(
      SELECT
            m.old_db_name,
            m.old_release,
            m.old_assembly,
            sie.old_version,
            sie.score
      FROM  stable_id_event sie, mapping_session m
      WHERE sie.mapping_session_id = m.mapping_session_id
      AND   sie.old_stable_id = "@{[$arch_id->stable_id]}"
      AND   m.old_db_name = "@{[$arch_id->db_name]}"
      $EXTRA_SQL
    );
  }

  my $sth = $self->prepare($sql);
  $sth->execute();
  my ($db_name, $release, $assembly, $version, $score) = $sth->fetchrow_array();
  $sth->finish();
  
  # you might have missed a stable ID that was deleted in the very first
  # stable ID mapping for this species, so go back and try again
  if( ! defined $db_name ) {
    my $sql_tmp = qq(
      SELECT
            m.old_db_name,
            m.old_release,
            m.old_assembly,
            sie.old_version,
            sie.score
      FROM  stable_id_event sie, mapping_session m
      WHERE sie.mapping_session_id = m.mapping_session_id
      AND   old_stable_id = "@{[$arch_id->stable_id]}"
      AND   m.old_db_name <> 'ALL'
      $EXTRA_SQL
      ORDER BY m.created DESC
    );
    $sql = $self->dbc->add_limit_clause($sql_tmp, 1);

    $sth = $self->prepare($sql);
    $sth->execute();
    ($db_name, $release, $assembly, $version, $score) = $sth->fetchrow_array();
    $sth->finish();
  }
  
  if( ! defined $db_name ) {
    # couldn't find stable ID in archive
    return 0;
  } else {
    $arch_id->version($version);
    if( ! defined $arch_id->{'db_name'} ) {
      $arch_id->db_name($db_name);
      $arch_id->release($release);
      $arch_id->assembly($assembly);
      $arch_id->score($score);
    }
  }

  return 1;
}


# static helper
sub _resolve_type {
  my $arch_id = shift;
  my $stable_id = $arch_id->stable_id();
  my $id_type;

  if( $stable_id =~ /.*G\d+$/ ) {
    $id_type = "Gene";
  } elsif( $stable_id =~ /.*T\d+$/ ) { 
    $id_type = "Transcript";
  } elsif( $stable_id =~ /.*P\d+$/ ) {
    $id_type = "Translation";
  } elsif( $stable_id =~ /.*E\d+$/ ) { 
    $id_type = "Exon";
  } else {
    $id_type = undef;
  }
  
  $arch_id->type( $id_type );
}


# deprecated methods (changed to more descriptive names)

sub fetch_pre_by_arch_id {
  my $self = shift;
  my $arch_id = shift;

  deprecate("Use fetch_predecessors_by_archive_id() instead");

  return $self->fetch_predecessors_by_archive_id($arch_id);
}

sub fetch_succ_by_arch_id {
  my $self = shift;
  my $arch_id = shift;

  deprecate("Use fetch_successors_by_archive_id() instead");

  return $self->fetch_successors_by_archive_id($arch_id);
}

sub fetch_all_currently_related {
  my $self = shift;
  my $arch_id = shift;

  deprecate("Use fetch_successor_history() instead");

  return $self->fetch_successor_history($arch_id);
}


=head2 fetch_all_by_gene_archive_id

  Description : DEPRECATED. Please use fetch_all_by_archive_id($archive_id,
                'Transcript') instead.

=cut

sub fetch_all_by_gene_archive_id {
  my $self = shift;
  my $gene_archive_id = shift;

  deprecate("Please use fetch_all_by_archive_id(\$archive_id, 'Transcript') instead.");

  return $self->fetch_all_by_archive_id($gene_archive_id, 'Transcript');
}

=head2 fetch_by_transcript_archive_id

  Description : DEPRECATED. Please use fetch_all_by_archive_id($archive_id,
                'Translation') instead. Note that the return type of the new
                method is different from this one.

=cut

sub fetch_by_transcript_archive_id {
  my $self = shift;
  my $transcript_archive_id = shift;

  deprecate("Please use fetch_all_by_archive_id(\$archive_id, 'Translation') instead.");

  my ($translation_archive_id) =
    @{ $self->fetch_all_by_archive_id($transcript_archive_id, 'Translation') };

  if ($translation_archive_id) {
    return $translation_archive_id;
  } else {
    return undef;
  }
}



1;


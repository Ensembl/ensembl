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
    fetch_all_by_gene_archive_id
    fetch_by_transcript_archive_id
    fetch_predecessors_by_archive_id
    fetch_successors_by_archive_id
    fetch_stable_id_history
    fetch_predecessor_history
    fetch_successor_history
    get_peptide
    list_dbnames
    _lookup_version
    _resolve_type

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
use Bio::EnsEMBL::Utils::Exception qw(deprecate);


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
  Description : Create an archiveStableId with given version and stableId
                No lookup is done in the database.
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


=head2 fetch_all_by_gene_archive_id

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId $gene_archive_id
  Example     : none
  Description : Given the ArchiveStableId of a gene retrieves ArchiveStableIds
                of Transcripts that make that gene.
  Returntype  : listref Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : Bio::EnsEMBL::ArchiveStableId->get_all_transcript_archive_ids
  Status      : At Risk
              : under development

=cut

sub fetch_all_by_gene_archive_id {
  my $self = shift;
  my $gene_archive_id = shift;
  my @result = ();

  my $sql = qq(
    SELECT
          ga.transcript_stable_id,
          ga.transcript_version,
          m.old_db_name,
          m.old_release,
          m.old_assembly
    FROM  gene_archive ga, mapping_session m
    WHERE ga.gene_stable_id = ?
    AND   ga.gene_version = ?
    AND   ga.mapping_session_id = m.mapping_session_id
  );
  
  my $sth = $self->prepare($sql);
  $sth->bind_param(1,$gene_archive_id->stable_id,SQL_VARCHAR);
  $sth->bind_param(2,$gene_archive_id->version,SQL_SMALLINT);
  $sth->execute();
  
  my ($stable_id, $version, $db_name, $release, $assembly);
  $sth->bind_columns(\$stable_id, \$version, \$db_name, \$release, \$assembly);

  while( $sth->fetch() ) {
    my $new_arch_id = Bio::EnsEMBL::ArchiveStableId->new
      (
       -version => $version,
       -adaptor => $self,
       -stable_id => $stable_id,
       -type => "Transcript",
       -db_name => $db_name,
       -release => $release,
       -assembly => $assembly
      );

    push( @result, $new_arch_id );
  }

  $sth->finish();
  return \@result;
}


=head2 fetch_by_transcript_archive_id

  Arg [1]     : Bio::EnsEMBL::ArchiveStableId
  Example     : none
  Description : Given a Transcripts ArchiveStableId retrieves the
                Translations ArchiveStableId. 
  Returntype  : Bio::EnsEMBL::ArchiveStableId or undef if not in database
  Exceptions  : none
  Caller      : Bio::EnsEMBL::ArchiveStableId->get_all_translation_archive_ids
  Status      : At Risk
              : under development

=cut

sub fetch_by_transcript_archive_id {
  my $self = shift;
  my $transcript_archive_id = shift;

  my $sql = qq(
    SELECT
          ga.translation_stable_id,
          ga.translation_version,
          m.old_db_name,
          m.old_release,
          m.old_assembly
    FROM  gene_archive ga, mapping_session m
    WHERE ga.transcript_stable_id = ?
    AND   ga.transcript_version = ?
  );
  
  my $sth = $self->prepare( $sql );
  $sth->bind_param(1,$transcript_archive_id->stable_id,SQL_VARCHAR);
  $sth->bind_param(2,$transcript_archive_id->version,SQL_SMALLINT);
  $sth->execute();
  
  my ($stable_id, $version, $db_name, $release, $assembly) = $sth->fetchrow_array();
  
  $sth->finish();

  if( $db_name ) {
    my $new_arch_id = Bio::EnsEMBL::ArchiveStableId->new
      (
       -version => $version,
       -adaptor => $self,
       -stable_id => $stable_id,
       -type => "Translation",
       -db_name => $db_name,
       -release => $release,
       -assembly => $assembly
      );

    return $new_arch_id;
  } else {
    return undef;
  }
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
          m.old_db_name,
          m.old_release,
          m.old_assembly
    FROM  mapping_session m, stable_id_event sie
    WHERE sie.mapping_session_id = m.mapping_session_id
    AND   sie.new_stable_id = ?
    AND   m.new_db_name = ?	
  );

  my $sth = $self->prepare( $sql );
  $sth->bind_param(1,$arch_id->stable_id, SQL_VARCHAR);
  $sth->bind_param(2,$arch_id->db_name,SQL_VARCHAR);
  $sth->execute();
  
  my ($old_stable_id, $old_version, $old_db_name, $old_release, $old_assembly);
  $sth->bind_columns(\$old_stable_id, \$old_version, \$old_db_name, \$old_release, \$old_assembly);
  
  while( $sth->fetch() ) {
    if( defined $old_stable_id ) {
      my $old_arch_id = Bio::EnsEMBL::ArchiveStableId->new
	( 
	 -version => $old_version,
	 -stable_id => $old_stable_id,
	 -db_name => $old_db_name,
         -release => $old_release,
         -assembly => $old_assembly,
	 -adaptor => $self
	);
      _resolve_type( $old_arch_id );
      push( @result, $old_arch_id );
    }
  }
  $sth->finish();

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
  
  my ($new_stable_id, $new_version, $new_db_name, $new_release, $new_assembly);
  $sth->bind_columns(\$new_stable_id, \$new_version, \$new_db_name, \$new_release, \$new_assembly);
  
  while( $sth->fetch() ) {
    if( defined $new_stable_id ) {
      my $new_arch_id = Bio::EnsEMBL::ArchiveStableId->new
	( 
	 -version => $new_version,
	 -stable_id => $new_stable_id,
	 -db_name => $new_db_name,
         -release => $new_release,
         -assembly => $new_assembly,
	 -adaptor => $self
	);
        
      _resolve_type($new_arch_id);
      push( @result, $new_arch_id );
    }
  }
  $sth->finish();
  
  return \@result;
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
    @$old = @$new;
    push @result, @$new;
  }

  # filter duplicates
  my %unique = map { join(":", $_->stable_id, $_->version, $_->release) => $_ }
    @result;
  @result = values %unique;

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
    @$old = @$new;
    push @result, @$new;
  }

  # filter duplicates
  my %unique = map { join(":", $_->stable_id, $_->version, $_->release) => $_ }
    @result;
  @result = values %unique;

  return \@result;
}


=head2 list_dbnames

  Args        : none
  Example     : none
  Description : A list of available database names from the latest (current) to
                the oldest (ordered).
  Returntype  : listref string
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
    
    my $first = 1;

    $sth->bind_columns( \$old_db_name, \$new_db_name );
    while( $sth->fetch() ) {
      if( $old_db_name eq "ALL" ) {
	next;
      }
      if( $first ) {
	$self->{'dbnames'} = [];
	push( @{$self->{'dbnames'}}, $new_db_name );
	push( @{$self->{'dbnames'}}, $old_db_name );
	$first = 0;
      } else {
	push( @{$self->{'dbnames'}}, $new_db_name );
      }
    }
    $sth->finish();
  }

  return $self->{'dbnames'};
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
    " and sie.type = '@{[lc($arch_id->{'type'})]}'" : '';

  if( ! defined $arch_id->{'db_name'} ) {
    # latest version of this stable id
      my $sql_tmp = qq(
        SELECT
              m.new_db_name,
              m.new_release,
              m.new_assembly,
              sie.new_version
        FROM  stable_id_event sie, mapping_session m
        WHERE sie.mapping_session_id = m.mapping_session_id
        AND   new_stable_id = "@{[$arch_id->stable_id]}" $EXTRA_SQL
        ORDER BY m.created DESC
      );
      $sql = $self->dbc->add_limit_clause($sql_tmp, 1);
  } else {
    $sql = qq(
      SELECT
            m.old_db_name,
            m.old_release,
            m.old_assembly,
            sie.old_version
      FROM  stable_id_event sie, mapping_session m
      WHERE sie.mapping_session_id = m.mapping_session_id
      AND   sie.old_stable_id = "@{[$arch_id->stable_id]}"
      AND   m.old_db_name = "@{[$arch_id->db_name]}"
      $EXTRA_SQL
    );
  }

  my $id_type;

  my $sth = $self->prepare( $sql );
  $sth->execute();
  my ($db_name, $release, $assembly, $version) = $sth->fetchrow_array();
  $sth->finish();
  
  if( ! defined $db_name ) {
    return 0;
  } else {
    $arch_id->version($version);
    if( ! defined $arch_id->{'db_name'} ) {
      $arch_id->db_name($db_name);
      $arch_id->release($release);
      $arch_id->assembly($assembly);
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


1;

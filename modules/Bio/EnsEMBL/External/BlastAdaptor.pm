=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::External::BlastAdaptor;

use strict;
use DBI;
use Storable qw(freeze thaw);
use Data::Dumper qw( Dumper );
use Time::Local;

use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::Search::HSP::EnsemblHSP; # This is a web module

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );
#@ISA = qw( Bio::EnsEMBL::DBSQL::DBAdaptor );


#----------------------------------------------------------------------
# Define SQL

#--- CREATE TABLES ---
our $SQL_CREATE_TICKET = "
CREATE TABLE blast_ticket (
  ticket_id int(10) unsigned NOT NULL auto_increment,
  create_time datetime NOT NULL default '0000-00-00 00:00:00',
  update_time datetime NOT NULL default '0000-00-00 00:00:00',
  ticket varchar(32) NOT NULL default '',
  status enum('CURRENT','DELETED') NOT NULL default 'CURRENT',
  object longblob,
  PRIMARY KEY  (ticket_id),
  UNIQUE KEY ticket (ticket),
  KEY create_time (create_time),
  KEY update_time (update_time)
) ENGINE=InnoDB";

our $SQL_CREATE_TABLE_LOG = "
CREATE TABLE blast_table_log (
  table_id int(10) unsigned NOT NULL auto_increment,
  table_name varchar(32),
  table_type enum('TICKET','RESULT','HIT','HSP') default NULL,
  table_status enum('CURRENT','FILLED','DELETED') default NULL,
  use_date date default NULL,
  create_time datetime default NULL,
  delete_time datetime default NULL,
  num_objects int(10) default NULL,
  PRIMARY KEY  (table_id),
  KEY table_name (table_name),
  KEY table_type (table_type),
  KEY use_date (use_date),
  KEY table_status (table_status)
) ENGINE=InnoDB";


our $SQL_CREATE_DAILY_RESULT = "
CREATE TABLE %s (
  result_id int(10) unsigned NOT NULL auto_increment,
  ticket varchar(32) default NULL,
  object longblob,
  PRIMARY KEY  (result_id),
  KEY ticket (ticket)
) ENGINE=InnoDB";

our $SQL_CREATE_DAILY_HIT = "
CREATE TABLE %s (
  hit_id int(10) unsigned NOT NULL auto_increment,
  ticket varchar(32) default NULL,
  object longblob,
  PRIMARY KEY  (hit_id),
  KEY ticket (ticket)
) ENGINE=InnoDB";

our $SQL_CREATE_DAILY_HSP = "
CREATE TABLE %s (
  hsp_id int(10) unsigned NOT NULL auto_increment,
  ticket varchar(32) default NULL,
  object longblob,
  chr_name varchar(32) default NULL,
  chr_start int(10) unsigned default NULL,
  chr_end int(10) unsigned default NULL,
  PRIMARY KEY  (hsp_id),
  KEY ticket (ticket)
) ENGINE=InnoDB MAX_ROWS=705032704 AVG_ROW_LENGTH=4000";

#--- TABLE LOG ---
our $SQL_SELECT_TABLE_LOG_CURRENT = "
SELECT   use_date
FROM     blast_table_log
WHERE    table_type   = ?
AND      table_status = 'CURRENT'
ORDER BY use_date DESC";

our $SQL_TABLE_LOG_INSERT = "
INSERT into blast_table_log 
       ( table_name, table_status, table_type, use_date, create_time)
VALUES ( ?, ?, ?, ?, NOW() )";

our $SQL_TABLE_LOG_UPDATE = "
UPDATE blast_table_log
SET    table_status = ?,
       delete_time  = ?,
       num_objects  = ?
WHERE  table_name   = ?";

#--- TICKETS ---

our $SQL_SEARCH_MULTI_STORE = "
INSERT INTO blast_ticket ( create_time, update_time, object, ticket )
VALUES                   ( NOW(), NOW(), ? , ? )";

our $SQL_SEARCH_MULTI_UPDATE = "
UPDATE blast_ticket
SET    object      = ?,
       update_time = NOW()
WHERE  ticket      = ?";

our $SQL_SEARCH_MULTI_RETRIEVE = "
SELECT object
FROM   blast_ticket
WHERE  ticket = ? ";

#--- RESULTS ---

our $SQL_RESULT_STORE = "
INSERT INTO blast_result%s ( object, ticket )
VALUES                   ( ? , ? )";

our $SQL_RESULT_UPDATE = "
UPDATE  blast_result%s
SET     object = ?,
        ticket = ?
WHERE   result_id = ?";

our $SQL_RESULT_RETRIEVE = "
SELECT object
FROM   blast_result%s
WHERE  result_id = ? ";

our $SQL_RESULT_RETRIEVE_TICKET = "
SELECT object
FROM   blast_result%s
WHERE  ticket = ? ";

#--- HITS ---

our $SQL_HIT_STORE = "
INSERT INTO blast_hit%s ( object, ticket )
VALUES                  ( ? , ? )";

our $SQL_HIT_UPDATE = "
UPDATE  blast_hit%s
SET     object = ?,
        ticket = ?
WHERE   hit_id = ?";

our $SQL_HIT_RETRIEVE = "
SELECT object
FROM   blast_hit%s
WHERE  hit_id = ? ";

#--- HSPS ---

our $SQL_HSP_STORE = "
INSERT INTO blast_hsp%s ( object, ticket, chr_name, chr_start, chr_end )
VALUES                  ( ? , ? , ? , ? , ? )";

our $SQL_HSP_UPDATE = "
UPDATE  blast_hsp%s
SET     object    = ?,
        ticket    = ?,
        chr_name  = ?,
        chr_start = ?,
        chr_end   = ?
WHERE   hsp_id    = ?";

our $SQL_HSP_RETRIEVE = "
SELECT object
FROM   blast_hsp%s
WHERE  hsp_id = ? ";

our $SQL_HSP_REMOVE = "
UPDATE  blast_hsp%s
SET     chr_name  = NULL,
        chr_start = NULL,
        chr_end   = NULL
WHERE   hsp_id    = ?";


#=head2 new
# 
#  Arg [1]   :
#  Function  :
#  Returntype:
#  Exceptions:
#  Caller    :
#  Example   :
# 
#=cut
#                                                                           
#
sub new {
  my $caller = shift;
#warn "DB - @_";
  my $connection = Bio::EnsEMBL::DBSQL::DBConnection->new(@_);
  my $self = $caller->SUPER::new($connection);
  $self->{'disconnect_flag'} = 1;
  return $self;
}
 

sub new_fast{
  my ($caller,$connection) = @_;
  my $self = $caller->SUPER::new($connection);
  $self->{'disconnect_flag'} = 1;
  return $self;
}

#----------------------------------------------------------------------

sub species {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_species} = $arg );
  $self->{_species};
}

#----------------------------------------------------------------------

=head2 ticket

  Arg [1]   : string ticket (optional)
  Function  : Get/get the blast ticket attribute
  Returntype: string ticket
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub ticket{
  my $key = "_ticket";
  my $self = shift;
  if( @_ ){ $self->{$key} = shift }
  return $self->{$key};
}

#----------------------------------------------------------------------

=head2 store

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub store {
  my $self = shift;
  my $obj = shift;
  my $ret_value = undef;
  if( $obj->isa("Bio::Tools::Run::SearchMulti") ) {
    $ret_value = $self->store_search_multi( $obj, @_ );
#    warn "Just stored as Bio::Tools::Run::SearchMulti";
  } elsif( $obj->isa( "Bio::Search::Result::ResultI" ) ) {
    $ret_value = $self->store_result(       $obj, @_ );
#    warn "Just stored as Bio::Tools::Result::ResultI";
  } elsif( $obj->isa( "Bio::Search::Hit::HitI" ) ) {
    $ret_value = $self->store_hit(       $obj, @_ );
#    warn "Just stored as Bio::Tools::Hit::HitI";
  } elsif( $obj->isa( "Bio::Search::HSP::HSPI" ) ) {
    $ret_value = $self->store_hsp(       $obj, @_ );
#    warn "Just stored as Bio::Tools::HSP::HSPI";
  } else {
#    warn "DID NOT STORE  ".ref($obj);
    $self->throw( "Do not know how to store objects of type ".ref($obj) );
    return undef;
  }
#  if( $self->{'disconnect_flag'} ) {
#    warn "HERE WE ARE DISCONNECTING....";
#    $self->dbc->db_handle->disconnect();
#    $self->dbc->connected(0);
#    warn "AND  WE ARE RECONNECTING....";
#    $self->dbc->connect();
#  }
  return $ret_value;
}

sub prepare {
  my $self = shift;
#  warn( "==> ", $self->dbc->dbname, " ", $self->dbc->db_handle );
#warn @_;
  my $T = $self->SUPER::prepare( @_ );
#  warn( "<== ", $self->dbc->dbname, " ", $self->dbc->db_handle );
  return $T;
}
#----------------------------------------------------------------------

=head2 retrieve

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub retrieve {
  my $self = shift;
  my $caller = shift;
  my %METHODS = qw(
    Bio::Tools::Run::EnsemblSearchMulti search_multi
    Bio::Search::Result::ResultI        result
    Bio::Search::Hit::HitI              hit
    Bio::Search::HSP::HSPI              hsp
  );
  foreach my $type (keys %METHODS) {
    if( UNIVERSAL::isa($caller, $type) ) {
      my $method = "retrieve_$METHODS{$type}";
      return $self->$method( @_ );
    }
  }
  return undef if UNIVERSAL::isa($caller,'Bio::Tools::Run::Search');
  $self->throw( "Do not know how to retrieve objects of type ". 
                ( ref($caller)? ref($caller) : $caller ) );
}

#----------------------------------------------------------------------

=head2 remove

  Arg [1]   : 
  Function  : TODO: implement remove functions
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub remove {
  my $self = shift;
  my $obj = shift;
  return 1 if $obj->isa("Bio::Tools::Run::EnsemblSearchMulti"); # Nothing to do here { return $self->remove_search_multi( @_ ); }
  return 1 if $obj->isa("Bio::Search::Result::ResultI");        # Nothing to do here { return $self->remove_result( @_ ); }
  return 1 if $obj->isa("Bio::Search::Hit::HitI");              # Nothing to do here { return $self->remove_hit( @_ ); }
  return $self->remove_hsp( $obj ) if $obj->isa("Bio::Search::HSP::HSPI");
  return undef(); # Do not know how to remove objects of this type 
}
#----------------------------------------------------------------------
=head2 store_search_multi

  Arg [1]   : Bio::Tools::Run::EnsemblSearchMulti obj
  Function  : Stores the ensembl SearchMulti container object in the database
  Returntype: scalar (token)
  Exceptions: 
  Caller    : 
  Example   : my $container_token = $blast_adpt->store_ticket( $container );

=cut

sub store_search_multi{
  my $self         = shift;
  my $search_multi = shift || 
    $self->throw( "Need a Bio::Tools::Run::EnsemblSearchMulti obj" );

  my $frozen = shift || $search_multi->serialise;

  my $dbh  = $self->dbc->db_handle;
  my $ticket  = $search_multi->token || $self->throw( "Bio::Tools::Run::EnsemblSearchMulti obj has no ticket" );

  my $sth = $self->prepare( $SQL_SEARCH_MULTI_RETRIEVE );
  my $rv = $sth->execute( $ticket ) ||  $self->throw( $sth->errstr );
  $sth->finish;

  if( $rv < 1 ){ # Insert (do first to minimise risk of race)
    my $sth = $self->prepare( $SQL_SEARCH_MULTI_STORE );
    $sth->execute( $frozen, $ticket ) || $self->throw( $sth->errstr );
    #$search_multi->token( $self->dbh->{mysql_insertid} );
    $sth->finish;
  }
  else{ # Update
    my $sth = $self->prepare( $SQL_SEARCH_MULTI_UPDATE );
    $sth->execute( $frozen, $ticket ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  my $sth = $self->prepare('show tables'); $sth->execute(); $sth->finish;
  return $search_multi->token();
}

#----------------------------------------------------------------------

=head2 retrieve_search_multi

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub retrieve_search_multi {
  my $self   = shift;
  my $ticket = shift || $self->throw( "Need an EnsemblSearchMulti ticket" );  

  my $dbh  = $self->dbc->db_handle;
warn $dbh;
warn $SQL_SEARCH_MULTI_RETRIEVE;
  my $sth = $self->prepare( $SQL_SEARCH_MULTI_RETRIEVE );
warn $sth;
  my $rv  = $sth->execute( $ticket ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $ticket not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
  $frozen || $self->throw( "Object from ticket $ticket is empty" );
  $sth->finish;
  return $frozen;
}



#----------------------------------------------------------------------
=head2 store_result

  Arg [1]   : Bio::Search::Result::EnsemblResult obj
  Function  : Stores the ensembl Result in the database
  Returntype: scalar (token)
  Exceptions: 
  Caller    : 
  Example   : my $result_token = $blast_adpt->store_result( $result );

=cut

sub store_result{
  my $self   = shift;
  my $res    = shift || $self->throw( "Need a Bio::Search::Result::EnsemblResult obj" );
  my $frozen = shift || $res->serialise;
  my $dbh    = $self->dbc->db_handle;
  my $sth;

  my ( $id, $use_date ) = split( '!!', $res->token || '' );
  $use_date ||= $self->use_date( 'RESULT' );
          #my $ticket = $res->group_ticket || warn( "Result $id has no ticket" );
  my $ticket  = $self->ticket || warn("Result $id BlastAdaptor has no ticket");

  my $rv = 0;
  if( $id ){
    $sth = $self->prepare( sprintf $SQL_RESULT_RETRIEVE, $use_date );
    $rv = $sth->execute( $id ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }
  if( $rv < 1 ){ # We have no result with this token string Insert
    my $use_date = $res->use_date() || $res->use_date($self->use_date('RESULT'));
    $sth = $self->prepare( sprintf $SQL_RESULT_STORE, $use_date );
    $sth->execute( $frozen, $ticket ) || $self->throw( $sth->errstr );
    my $id = $dbh->{mysql_insertid};
    $res->token( join( '!!', $id, $use_date ) );
    $sth->finish;
  } else {  # Update
    $sth = $self->prepare( sprintf $SQL_RESULT_UPDATE, $use_date );
    $sth->execute( $frozen, $ticket, $id ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  return $res->token();
}

sub store_result_2{
  my $self   = shift;
  my $res    = shift || $self->throw( "Need a Bio::Search::Result::EnsemblResult obj" );
  my $frozen = shift || $res->serialise;
  my $dbh    = $self->dbc->db_handle;
  my $sth;

  my ( $id, $use_date ) = split( '!!', $res->token || '' );
  $use_date ||= $self->use_date( 'RESULT' );
          #my $ticket = $res->group_ticket || warn( "Result $id has no ticket" );
  my $ticket  = $self->ticket || warn("Result $id BlastAdaptor has no ticket");

  my $rv = 0;
  if( $ticket ){
    $sth = $self->prepare( sprintf $SQL_RESULT_RETRIEVE_TICKET, $use_date );
    $rv = $sth->execute( $ticket ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  } 
  if( !$rv && $id ){
    $sth = $self->prepare( sprintf $SQL_RESULT_RETRIEVE, $use_date );
    $rv = $sth->execute( $id ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }
  if( $rv < 1 ){ # We have no result with this token string Insert
    my $use_date = $res->use_date() || $res->use_date($self->use_date('RESULT'));
    $sth = $self->prepare( sprintf $SQL_RESULT_STORE, $use_date );
    $sth->execute( $frozen, $ticket ) || $self->throw( $sth->errstr );
    my $id = $dbh->{mysql_insertid};
    $res->token( join( '!!', $id, $use_date ) );
    $sth->finish;
  } else {  # Update
    $sth = $self->prepare( sprintf $SQL_RESULT_UPDATE, $use_date );
    $sth->execute( $frozen, $ticket, $id ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  return $res->token();
}

#----------------------------------------------------------------------

=head2 retrieve_result

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub retrieve_result{
  my $self = shift;
  my $token  = shift || $self->throw( "Need a Hit token" );
  my ( $id, $use_date ) = split( '!!',$token);
  $use_date ||= '';

  my $dbh  = $self->dbc->db_handle;
  my $sth = $self->prepare( sprintf $SQL_RESULT_RETRIEVE, $use_date );
  my $rv  = $sth->execute( $id ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $id not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
  $frozen || $self->throw( "Object from result $id is empty" );
  $sth->finish;
  return $frozen;
}

#----------------------------------------------------------------------
=head2 store_hit

  Arg [1]   : Bio::Search::Hit::EnsemblHit obj
  Function  : Stores the ensembl Hit in the database
  Returntype: scalar (token)
  Exceptions: 
  Caller    : 
  Example   : my $hit_token = $blast_adpt->store_hit( $hit );

=cut

sub store_hit{
  my $self = shift;
  my $hit  = shift || 
    $self->throw( "Need a Bio::Search::Hit::EnsemblHit obj" );
  my $frozen = shift || $hit->serialise;

  my $dbh  = $self->dbc->db_handle;

  my ( $id, $use_date ) = split( '!!', $hit->token || '' );
  $use_date ||= $hit->use_date() || $hit->use_date($self->use_date('HIT'));;
  #my $ticket = $hit->group_ticket || warn( "Hit $id has no ticket" );
  my $ticket = $self->ticket || warn("Hit $id BlastAdaptor has no ticket");

  my $rv = 0;
  if( $id ){
    my $sth = $self->prepare( sprintf $SQL_HIT_RETRIEVE, $use_date );
    $rv = $sth->execute( $id ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }
  if( $rv < 1 ){ # Insert
    my $sth = $self->prepare( sprintf $SQL_HIT_STORE, $use_date );
    $sth->execute( $frozen, $ticket ) || $self->throw( $sth->errstr );
    my $id = $dbh->{mysql_insertid};
    $hit->token( join( '!!', $id, $use_date ) );
    $sth->finish;
  }
  else{ # Update
    my $sth = $self->prepare( sprintf $SQL_HIT_UPDATE, $use_date );
    $sth->execute( $frozen, $ticket, $id ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  return $hit->token();
}

#----------------------------------------------------------------------

=head2 retrieve_hit

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub retrieve_hit{
  my $self   = shift;
  my $token  = shift || $self->throw( "Need a Hit token" );
  my ( $id, $use_date ) = split( '!!',$token);
  $use_date ||= '';
  my $dbh  = $self->dbc->db_handle;
  my $sth = $self->prepare( sprintf $SQL_HIT_RETRIEVE, $use_date );
  my $rv  = $sth->execute( $id ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $token not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
  $frozen || $self->throw( "Object from hit $id is empty" );
  $sth->finish;
  return $frozen;
}

#----------------------------------------------------------------------
=head2 store_hsp

  Arg [1]   : Bio::Search::HSP::EnsemblHSP obj
  Function  : Stores the ensembl HSP in the database
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub store_hsp{
  my $self = shift;
  my $hsp  = shift || 
    $self->throw( "Need a Bio::Search::HSP::EnsemblHSP obj" );
  my $frozen = shift || $hsp->serialise;

  my $dbh  = $self->dbc->db_handle;
  my ( $id, $use_date ) = split( '!!', $hsp->token || '');
  $use_date ||= $hsp->use_date() || $hsp->use_date($self->use_date('HSP'));

  #my $ticket = $hsp->group_ticket || warn( "HSP $id has no ticket" );
  my $ticket = $self->ticket || warn( "HSP $id BlastAdaptor has no ticket" );

  my $chr_name  = '';
  my $chr_start = 0;
  my $chr_end   = 0;
  if( my $genomic = $hsp->genomic_hit ){
    $chr_name  = $genomic->seq_region_name;
    $chr_start = $genomic->start;
    $chr_end   = $genomic->end;
  }
  my $rv = 0;
  if( $id ){
    my $sth = $self->prepare( sprintf $SQL_HSP_RETRIEVE, $use_date );
    $rv = $sth->execute( $id ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }
  if( $rv < 1 ){ # Insert
    my $use_date = $hsp->use_date() || $hsp->use_date($self->use_date('HSP'));
    my $sth = $self->prepare( 'show tables' ); $sth->execute(); $sth->finish();
    $sth = $self->prepare( sprintf $SQL_HSP_STORE, $use_date );
    my @bound = ( $frozen, $ticket, $chr_name,  $chr_start, $chr_end );
    $sth->execute( @bound ) || $self->throw( $sth->errstr );
    my $id = $dbh->{mysql_insertid};
    $hsp->token( join( '!!', $id, $use_date ) );
    $sth->finish;
  }
  else{ # Update
    my $sth = $self->prepare( sprintf $SQL_HSP_UPDATE, $use_date );
    my @bound = ( $frozen, $ticket, $chr_name,  $chr_start, $chr_end, $id );
    $sth->execute( @bound ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  return $hsp->token();
}

#----------------------------------------------------------------------

=head2 retrieve_hsp

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub retrieve_hsp{
  my $self   = shift;
  my $token  = shift || $self->throw( "Need an HSP token" );
  my ( $id, $use_date ) = split( '!!',$token);
  $use_date ||= '';
  my $dbh  = $self->dbc->db_handle;
  my $sth = $self->prepare( sprintf $SQL_HSP_RETRIEVE, $use_date );
  my $rv  = $sth->execute( $id ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $token not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
  $frozen || $self->throw( "Object from hsp $id is empty" );
  $sth->finish;
  return $frozen;
}

#----------------------------------------------------------------------

=head2 remove_hsp

  Arg [1]   : $hsp object to be removed
  Function  : 'removes' hsp from e.g. contigview by setting chr fields
              to null
  Returntype: 
  Exceptions: 
  Caller    : $self->remove
  Example   : 

=cut

sub remove_hsp {
  my $self = shift;
  my $hsp  = shift || 
    $self->throw( "Need a Bio::Search::HSP::EnsemblHSP obj" );

  my $dbh  = $self->dbc->db_handle;

  my ( $id, $use_date ) = split( '!!', $hsp->token || '');
  $use_date ||= $hsp->use_date() || $hsp->use_date($self->use_date('HSP'));

  my $sth = $self->prepare( sprintf $SQL_HSP_REMOVE, $use_date );
  my @bound = ( $id );
  my $rv = $sth->execute( @bound ) || $self->throw( $sth->errstr );
  $sth->finish;
  return 1;
}



#----------------------------------------------------------------------

=head2 get_all_HSPs

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub get_all_HSPs {
   my $self = shift;
   my $ticket    = shift || $self->throw( "Need a search ticket!");
   my $chr_name  = shift || undef;
   my $chr_start = shift || undef;
   my $chr_end   = shift || undef;
   my ( $id, $use_date )  = split( '!!', $ticket );
   $use_date ||= '';

   my $SQL = qq(
SELECT object, hsp_id
FROM   blast_hsp%s
WHERE  ticket = ? );

   my $CHR_SQL = qq(
AND    chr_name = ? );

   my $RANGE_SQL = qq(
AND    chr_start <= ?
AND    chr_end   >= ? );

   my $q = sprintf( $SQL, $use_date );
   my @binded = ( $id );

   if( $chr_name ){
     $q .= $CHR_SQL;
     push @binded, $chr_name;

     if( $chr_start && $chr_end ){
       $q .= $RANGE_SQL;
       push @binded, $chr_end, $chr_start;
     }
   }
   my $sth = $self->dbc->db_handle->prepare($q);
   my $rv = $sth->execute( @binded ) || $self->throw( $sth->errstr );

   my @hsps = ();
   foreach my $row( @{$sth->fetchall_arrayref()} ){
     # Retrieve HSP and reset token
     my $hsp = thaw( $row->[0] );
     my $hsp_id = $row->[1];
     $hsp->token( join( '!!', $hsp_id, $use_date  ) );
     push @hsps, $hsp;
   }
   $sth->finish;
   return [@hsps];
}



#----------------------------------------------------------------------

=head2 get_all_SearchFeatures

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub get_all_SearchFeatures {
  my $self = shift;
  my $hsps = $self->get_all_HSPs(@_);
  my $ticket = shift;

  $self->dynamic_use( ref($hsps->[0] ) );
  my @feats = ();
  foreach my $hsp( @$hsps ){
    my $base_align = $hsp->genomic_hit || next;

    ( $ticket ) = split( "!!", $ticket );
    my $hsp_id = join( "!!", $ticket, $hsp->token );

    $base_align->hseqname( join( ":", $base_align->hseqname, $hsp_id ) );
    push @feats, $base_align;
  }
  return [ @feats ];
}

sub dynamic_use {
  my( $self, $classname ) = @_;
  my( $parent_namespace, $module ) = $classname =~/^(.*::)(.*?)$/;
  no strict 'refs';
  return 1 if $parent_namespace->{$module.'::'}; # return if already used
  eval "require $classname";
  if($@) {
    warn "DrawableContainer: failed to use $classname\nDrawableContainer: $@";
    return 0;
  }
  $classname->import();
  return 1;
}

#----------------------------------------------------------------------

=head2 use_date

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

my %valid_table_types = ( HIT=>1, HSP=>1, RESULT=>1 );
sub use_date {
  my $key  = '_current_table';
  my $self = shift;
  my $type = uc( shift );
#warn "$self --- $key --- $type $self";
  $valid_table_types{$type} || 
    $self->throw( "Need a table type (Result, Hit or HSP)" );

  $self->{$key} ||= {};
  if( ! $self->{$key}->{$type} ){
    my $sth = $self->dbc->db_handle->prepare( "
SELECT table_type,  use_date
  FROM blast_table_log
 WHERE table_status = 'CURRENT'
ORDER BY use_date ASC" );
#warn "prepare... $sth";
#warn $SQL_SELECT_TABLE_LOG_CURRENT;
#warn $type;
    my $rv = $sth->execute();# $type );
#warn $rv;
    unless( $rv ) { 
      $sth->finish;
      warn( $sth->errstr );
      return;
    }
#warn "exec...";
    foreach my $r (@{ $sth->fetchall_arrayref }) {
      my $date = $r->[1];
      $date =~ s/-//g;
      $self->{$key}->{$r->[0]} = $date;
#warn "$r->[0] ---> $r->[1] ---> $date";
    }
#    $rv > 0 || ( warn( "No current $type table found" ) && return );
    $sth->finish;
#warn "end of finish...";
  }
  return $self->{$key}->{$type};
}



#----------------------------------------------------------------------

=head2 clean_blast_database

  Arg [1]   : int $days
  Function  : Removes blast tickets older than $days days
  Returntype: 
  Exceptions: SQL errors
  Caller    : 
  Example   : $ba->clean_blast_database(14)

=cut

sub clean_blast_database{
  my $self = shift;
  my $days = shift || $self->throw( "Missing arg: number of days" );
  $days =~ /\D/    && $self->throw( "Bad arg: number of days $days not int" );
  my $dbh = $self->dbc->db_handle;

  # Get list of tickets > $days days old
  my $q = qq(
    SELECT ticket_id
      FROM blast_ticket
     WHERE update_time < SUBDATE( NOW(), INTERVAL $days DAY ) );

  my $sth = $self->dbc->db_handle->prepare($q);
  my $rv = $sth->execute() || $self->throw( $sth->errstr );
  my $res = $sth->fetchall_arrayref;
  $sth->finish;

  # Delete result and ticket rows associated with old tickets
  my $q_del_tmpl = qq(
     DELETE
       FROM blast_ticket
      WHERE ticket_id = %s);

  my $c = 0;
  foreach my $row( @$res ){
    my $ticket_id = $row->[0];
    $c++;
    my $q_del = sprintf( $q_del_tmpl, $ticket_id );
    my $sth = $self->dbc->db_handle->prepare($q_del);
    my $rv = $sth->execute() || $self->throw( $sth->errstr );
  }
  warn "Purging $days days: Deleted $c rows\n";

  # Drop daily Result, Hit and HSP tables not updated within $days days
  my $q_find = 'show table status like ?';
  my $sth2 = $self->prepare( $q_find );
  $sth2->execute( "blast_result%" ) || $self->throw( $sth2->errstr );
  my $res_res = $sth2->fetchall_arrayref();
  $sth2->execute( "blast_hit%" ) || $self->throw( $sth2->errstr );
  my $hit_res = $sth2->fetchall_arrayref();
  $sth2->execute( "blast_hsp%" ) || $self->throw( $sth2->errstr );
  my $hsp_res = $sth2->fetchall_arrayref();

  my @deletable_hit_tables;
  foreach my $row( @$res_res, @$hit_res, @$hsp_res ){
    my $table_name  = $row->[0];  ## table name
    my $num_rows    = $row->[4];  ## # Rows...
    my $update_time = $row->[12]; ## update time ---  Should be a string like 2003-08-15 10:36:56

    #cope with an innodb tables that have no update time (keep them for five days longer)
    my $days_to_keep = $days;
    unless ($update_time) {
      $update_time = $row->[11];
      $days_to_keep = $days + 5;
      next unless $update_time;
    }


    my @time = split( /[-:\s]/, $update_time );
    my $epoch_then = timelocal( $time[5], $time[4],   $time[3], 
				$time[2], $time[1]-1, $time[0] - 1900 );
    my $secs_old = time() - $epoch_then;
    my $days_old = $secs_old / ( 60 * 60 * 24 );
    if( $days_old > $days_to_keep ){
      warn( "Dropping table $table_name: $num_rows rows\n" );
      my $sth_drop = $self->prepare( "DROP table $table_name" );
      my $sth_log  = $self->prepare( $SQL_TABLE_LOG_UPDATE );
      $sth_drop->execute || $self->throw( $sth_drop->errstr );
      my( $se,$mi,$hr,$da,$mo,$yr ) = (localtime)[0,1,2,3,4,5];
      my $now = sprintf( "%4d-%2d-%2d %2d:%2d:%2d", 
			 $yr+1900,$mo+1,$da,$hr,$mi,$se );
      $sth_log->execute
        ('DELETED',$now,$num_rows,$table_name) ||
	  $self->throw( $sth_log->errstr );
    }
  }

  return 1;
}

#----------------------------------------------------------------------

=head2 create_tables

  Arg [1]   : none
  Function  : Creates the blast_ticket and blast_table_log
              tables in the database indicated by the database handle.
              Checks first to make sure they do not exist
  Returntype: boolean
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub create_tables {
  my $self = shift;
  my $dbh = $self->dbc->db_handle;

  # Get list of existing tables in database
  my $q = 'show tables like ?';
  my $sth = $self->prepare( $q );
  my $rv_tck = $sth->execute("blast_ticket")    || $self->throw($sth->errstr);
  my $rv_log = $sth->execute("blast_table_log" )|| $self->throw($sth->errstr);
  $sth->finish;

  if( $rv_tck == 0 ){
    warn( "Creating blast_ticket table\n" );
    my $sth = $self->prepare( $SQL_CREATE_TICKET );
    my $rv = $sth->execute() || $self->throw( $sth->errstr );
    $sth->finish;
  }
  else{ warn( "blast_ticket table already exists\n" ) }

  if( $rv_log == 0 ){
    warn( "Creating blast_result table\n" );
    my $sth = $self->prepare( $SQL_CREATE_TABLE_LOG );
    my $rv = $sth->execute() || $self->throw( $sth->errstr );    
     $sth->finish;
  }
  else{ warn( "blast_table_log table already exists\n" ) }  

  return 1;
}

#----------------------------------------------------------------------

=head2 rotate_daily_tables

  Arg [1]   : none
  Function  : Creates the daily blast_result{date}, blast_hit{date} 
              and blast_hsp{date} tables in the database indicated by 
              the database handle.
              Checks first to make sure they do not exist.
              Sets the new table to 'CURRENT' in the blast_table_log.
              Sets the previous 'CURRENT' table to filled.
  Returntype: boolean
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub rotate_daily_tables {
  my $self = shift;
  my $dbh = $self->dbc->db_handle;
 
  # Get date
  my( $day, $month, $year ) = (localtime)[3,4,5];
  my $date = sprintf( "%04d%02d%02d", $year+1900, $month+1, $day );

  my $res_table = "blast_result$date";
  my $hit_table = "blast_hit$date";
  my $hsp_table = "blast_hsp$date";

  # Get list of existing tables in database
  my $q = 'show table status like ?';
  my $sth = $self->prepare( $q );
  my $rv_res  = $sth->execute($res_table) || $self->throw($sth->errstr);
  my $rv_hit  = $sth->execute($hit_table) || $self->throw($sth->errstr);
  my $rv_hsp  = $sth->execute($hsp_table) || $self->throw($sth->errstr);
  $sth->finish;

  if( $rv_res == 0 ){
    warn( "Creating today's $res_table table\n" );

    # Create new table
    my $q = sprintf($SQL_CREATE_DAILY_RESULT, $res_table);
    my $sth1 = $self->prepare( $q );
    my $rv = $sth1->execute() || $self->throw( $sth1->errstr );

    # Flip current table in blast_table_tog
    my $last_date = $self->use_date( "RESULT" ) || '';
    my $sth2 = $self->prepare( $SQL_TABLE_LOG_INSERT );
    my $sth3 = $self->prepare( $SQL_TABLE_LOG_UPDATE );
    $sth2->execute( "$res_table",'CURRENT','RESULT',$date ) 
      || die( $self->throw( $sth2->errstr ) );
    $sth3->execute( 'FILLED','0',0,"blast_result$last_date") 
      || die( $self->throw( $sth3->errstr ) );
    $sth1->finish();
    $sth2->finish();
    $sth3->finish();    
  }
  else{ warn( "Today's $res_table table already exists\n" ) }

  if( $rv_hit == 0 ){
    warn( "Creating today's $hit_table table\n" );

    # Create new table
    my $q = sprintf($SQL_CREATE_DAILY_HIT, $hit_table);
    my $sth1 = $self->prepare( $q );
    my $rv = $sth1->execute() || $self->throw( $sth1->errstr );         

    # Flip current table in blast_table_tog
    my $last_date = $self->use_date( "HIT" ) || '';
    my $sth2 = $self->prepare( $SQL_TABLE_LOG_INSERT );
    my $sth3 = $self->prepare( $SQL_TABLE_LOG_UPDATE );
    $sth2->execute( "$hit_table",'CURRENT','HIT',$date ) 
      || die( $self->throw( $sth2->errstr ) );
    $sth3->execute( 'FILLED','0',0,"blast_hit$last_date") 
      || die( $self->throw( $sth3->errstr ) );
    $sth1->finish();
    $sth2->finish();
    $sth3->finish();    
  }
  else{ warn( "Today's $hit_table table already exists\n" ) }
  
  if( $rv_hsp == 0 ){
    warn( "Creating today's $hsp_table table\n" );

    # Create new table
    my $q = sprintf($SQL_CREATE_DAILY_HSP, $hsp_table );
    my $sth1 = $self->prepare( $q );
    my $rv = $sth1->execute() || $self->throw( $sth1->errstr );         

    # Flip current table in blast_table_tog
    my $last_date = $self->use_date( "HSP" ) || '';    
    my $sth2 = $self->prepare( $SQL_TABLE_LOG_INSERT );
    my $sth3 = $self->prepare(  $SQL_TABLE_LOG_UPDATE );
    $sth2->execute( "$hsp_table",'CURRENT','HSP',$date ) 
      || die( $self->throw( $sth2->errstr ) );
    $sth3->execute( 'FILLED','0',0,"blast_hsp$last_date") 
      || die( $self->throw( $sth3->errstr ) );
    $sth1->finish();
    $sth2->finish();
    $sth3->finish();    
  }
  else{ warn( "Today's $hsp_table table already exists\n" ) }
  return 1;
}

#----------------------------------------------------------------------


=head2 cleanup_processes

  Arg [1]   : none
  Function  : Kills any sleeping processes older that 1000
  Returntype: boolean
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub cleanup_processes {
  my $self = shift;
  my $dbh = $self->dbc->db_handle;
  my $sth = $self->prepare( 'show processlist' );
  my $kill_sth = $self->prepare('kill ?');
  $sth->execute;
  my $res = $sth->fetchall_arrayref([0,3,4,5]);
  my $c = 0;
  foreach my $ps (@$res) {
    my ($pid,$db,$stat,$time) = @$ps;
    if ($db eq 'ensembl_blast') {
      if ( ($stat eq 'Sleep') && ($time > 1000) ) {
	$kill_sth->execute($pid);
	$c++;
      }
    }
  }
  warn "Killed $c processes";
  return 1;
}




1;

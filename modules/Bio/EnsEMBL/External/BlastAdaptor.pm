# Let the code begin...
package Bio::EnsEMBL::External::BlastAdaptor;

use strict;
use DBI;
use Storable qw(freeze thaw);
use Data::Dumper qw( Dumper );

use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::Search::HSP::EnsemblHSP; # This is a web module

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );
#@ISA = qw( Bio::EnsEMBL::DBSQL::DBConnection );


#----------------------------------------------------------------------
# Define SQL

#--- CREATE TABLES ---
our $SQL_CREATE_TICKET = "
CREATE TABLE blast_ticket (
  ticket_id int(10) unsigned NOT NULL auto_increment,
  create_time datetime NOT NULL default '0000-00-00 00:00:00',
  update_time datetime NOT NULL default '0000-00-00 00:00:00',
  ticket varchar(32) NOT NULL default '',
  object longblob,
  PRIMARY KEY  (ticket_id),
  UNIQUE KEY ticket (ticket),
  KEY create_time (create_time),
  KEY update_time (update_time)
) TYPE=MyISAM";

our $SQL_CREATE_RESULT = "
CREATE TABLE blast_result (
  result_id int(10) unsigned NOT NULL auto_increment,
  ticket varchar(32) default NULL,
  object longblob,
  PRIMARY KEY  (result_id),
  KEY ticket (ticket)
) TYPE=MyISAM";

our $SQL_CREATE_TABLE_LOG = "
CREATE TABLE blast_table_log (
  table_id int(10) unsigned NOT NULL auto_increment,
  table_type enum('TICKET','RESULT','HIT','HSP') default NULL,
  table_status enum('CURRENT','FILLED','DELETED') default NULL,
  use_date date default NULL,
  create_time datetime default NULL,
  delete_time datetime default NULL,
  num_objects int(10) default NULL,
  PRIMARY KEY  (table_id),
  KEY table_type (table_type),
  KEY use_date (use_date),
  KEY table_status (table_status)
) TYPE=MyISAM";

our $SQL_CREATE_DAILY_HIT = "
CREATE TABLE blast_hit%s (
  hit_id int(10) unsigned NOT NULL auto_increment,
  ticket varchar(32) default NULL,
  object longblob,
  PRIMARY KEY  (hit_id),
  KEY ticket (ticket)
) TYPE=MyISAM";

our $SQL_CREATE_DAILY_HSP = "
CREATE TABLE blast_hsp%s (
  hsp_id int(10) unsigned NOT NULL auto_increment,
  ticket varchar(32) default NULL,
  object longblob,
  chr_name varchar(32) default NULL,
  chr_start int(10) unsigned default NULL,
  chr_end int(10) unsigned default NULL,
  PRIMARY KEY  (hsp_id),
  KEY ticket (ticket)
) TYPE=MyISAM MAX_ROWS=705032704 AVG_ROW_LENGTH=4000";

#--- TABLE LOG ---
our $SQL_SELECT_TABLE_LOG_CURRENT = "
SELECT   use_date
FROM     blast_table_log
WHERE    table_type   = ?
AND      table_status = 'CURRENT'
ORDER BY use_date DESC";

our $SQL_TABLE_LOG_INSERT = "
INSERT into blast_table_log 
       ( table_status, table_type, use_date, create_time)
VALUES ( ?, ?, ?, NOW() )";

our $SQL_TABLE_LOG_UPDATE = "
UPDATE blast_table_log
SET    table_status = ?
WHERE  table_type   = ?
AND    use_date     = ?";

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
INSERT INTO blast_result ( object, ticket )
VALUES                   ( ? , ? )";

our $SQL_RESULT_UPDATE = "
UPDATE  blast_result
SET     object = ?,
        ticket = ?
WHERE   result_id = ?";

our $SQL_RESULT_RETRIEVE = "
SELECT object
FROM   blast_result
WHERE  result_id = ? ";

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
        chr_name  = ?,
        chr_start = ?,
        chr_end   = ?
WHERE   hsp_id    = ?";

our $SQL_HSP_RETRIEVE = "
SELECT object
FROM   blast_hsp%s
WHERE  hsp_id = ? ";


#----------------------------------------------------------------------

=head2 new

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub new {
  my $caller = shift;
  my $connection = Bio::EnsEMBL::DBSQL::DBConnection->new(@_);
  my $self = $caller->SUPER::new($connection);
  return $self;
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

  my $dbh  = $self->db->db_handle;

  my $store_obj = $search_multi->_prepare_storable;
  my $frozen  = freeze( $store_obj );
  my $ticket  = $search_multi->token ||
    $self->throw( "Bio::Tools::Run::EnsemblSearchMulti obj has no ticket" );

  my $sth = $dbh->prepare( $SQL_SEARCH_MULTI_RETRIEVE );
  my $rv = $sth->execute( $ticket ) ||  $self->throw( $sth->errstr );
  $sth->finish;

  if( $rv < 1 ){ # Insert (do first to minimise risk of race)
    my $sth = $dbh->prepare( $SQL_SEARCH_MULTI_STORE );
    $sth->execute( $frozen, $ticket ) || $self->throw( $sth->errstr );
    #$search_multi->token( $self->dbh->{mysql_insertid} );
    $sth->finish;
  }
  else{ # Update
    my $sth = $dbh->prepare( $SQL_SEARCH_MULTI_UPDATE );
    $sth->execute( $frozen, $ticket ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
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

  my $dbh  = $self->db->db_handle;
  my $sth = $dbh->prepare( $SQL_SEARCH_MULTI_RETRIEVE );
  my $rv  = $sth->execute( $ticket ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $ticket not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
  $frozen || $self->throw( "Object from ticket $ticket is empty" );
  $sth->finish;

  my $stored_obj = thaw( $frozen );
  if( ! ref( $stored_obj ) or 
      ! $stored_obj->isa( 'Bio::Root::Storable' ) ){
    $self->throw( "Token $ticket returned no data" );
  }
  return $stored_obj;
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
  my $self = shift;
  my $res  = shift || 
    $self->throw( "Need a Bio::Search::Result::EnsemblResult obj" );

  my $dbh  = $self->db->db_handle;

  my $ticket = $res->group_ticket;
  my $token  = $res->token;
  my $store_obj = $res->_prepare_storable;
  my $frozen = freeze( $store_obj );

  my $rv = 0;
  if( $token ){
    my $sth = $dbh->prepare( $SQL_RESULT_RETRIEVE );
    $rv = $sth->execute( $token ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }
  if( $rv < 1 ){# Insert
    my $sth = $dbh->prepare( $SQL_RESULT_STORE );
    $sth->execute( $frozen, $ticket ) || $self->throw( $sth->errstr );
    $res->token( $dbh->{mysql_insertid} );
    $sth->finish;
  }
  else{  # Update
    my $sth = $dbh->prepare( $SQL_RESULT_UPDATE );
    $sth->execute( $frozen, $ticket, $token ) || $self->throw( $sth->errstr );
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
  my $id   = shift || $self->throw( "Need a Result token" );

  my $dbh  = $self->db->db_handle;
  my $sth = $dbh->prepare( $SQL_RESULT_RETRIEVE );
  my $rv  = $sth->execute( $id ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $id not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
  $frozen || $self->throw( "Object from result $id is empty" );
  $sth->finish;

  my $stored_obj = thaw( $frozen );
  if( ! ref( $stored_obj ) or 
      ! $stored_obj->isa( 'Bio::Root::Storable' ) ){
    $self->throw( "Token $id returned no data" );
  }
  return $stored_obj;
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

  my $dbh  = $self->db->db_handle;

  my $store_obj = $hit->_prepare_storable;
  my $frozen = freeze( $store_obj );
  my $ticket = $hit->group_ticket;
  my ( $id, $use_date ) = split( '!!', $hit->token);
  $use_date ||= '';

  my( $rv );
  if( $id ){
    my $sth = $dbh->prepare( sprintf $SQL_HIT_RETRIEVE, $use_date );
    $rv = $sth->execute( $id ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }
  if( $rv > 0 ){ # Update
    my $sth = $dbh->prepare( sprintf $SQL_HIT_UPDATE, $use_date );
    $sth->execute( $frozen, $ticket, $id ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  else{ # Insert
    my $use_date = $self->use_date('HIT') || '';
    my $sth = $dbh->prepare( sprintf $SQL_HIT_STORE, $use_date );
    $sth->execute( $frozen, $ticket ) || $self->throw( $sth->errstr );
    my $id = $dbh->{mysql_insertid};
    $hit->token( join( '!!', $id, $use_date ) );
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
  
  my $dbh  = $self->db->db_handle;
  my $sth = $dbh->prepare( sprintf $SQL_HIT_RETRIEVE, $use_date );
  my $rv  = $sth->execute( $id ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $token not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
  $frozen || $self->throw( "Object from hit $id is empty" );
  $sth->finish;

  my $stored_obj = thaw( $frozen );
  if( ! ref( $stored_obj ) or 
      ! $stored_obj->isa( 'Bio::Root::Storable' ) ){
    $self->throw( "Token $token returned no data" );
  }
  return $stored_obj;
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
  my $dbh  = $self->db->db_handle;

  my $store_obj = $hsp->_prepare_storable;
  my $frozen = freeze( $store_obj );
  my $ticket = $hsp->group_ticket;
  my ( $id, $use_date ) = split( '!!', $hsp->token);
  $use_date ||= '';

  my $chr_name  = 'NULL';
  my $chr_start = 'NULL';
  my $chr_end   = 'NULL';
  if( my $genomic = $hsp->genomic_hit ){
    $chr_name  = $genomic->seq_id;
    $chr_start = $genomic->start;
    $chr_end   = $genomic->end;
  } 

  my( $rv );
  if( $id ){
    my $sth = $dbh->prepare( sprintf $SQL_HSP_RETRIEVE, $use_date );
    $rv = $sth->execute( $id ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }
  if( $rv > 0 ){ # Update
    my $sth = $dbh->prepare( sprintf $SQL_HSP_UPDATE, $use_date );
    my @bound = ( $frozen, $ticket, $chr_name,  $chr_start, $chr_end, $id );
    $sth->execute( @bound ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  else{ # Insert
    my $use_date = $self->use_date('HSP') || '';
    my $sth = $dbh->prepare( sprintf $SQL_HSP_STORE, $use_date );
    my @bound = ( $frozen, $ticket, $chr_name,  $chr_start, $chr_end );
    $sth->execute( @bound ) || $self->throw( $sth->errstr );
    my $id = $dbh->{mysql_insertid};
    $hsp->token( join( '!!', $id, $use_date ) );
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

  my $dbh  = $self->db->db_handle;
  my $sth = $dbh->prepare( sprintf $SQL_HSP_RETRIEVE, $use_date );
  my $rv  = $sth->execute( $id ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $token not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
  $frozen || $self->throw( "Object from hsp $id is empty" );
  $sth->finish;

  my $stored_obj = thaw( $frozen );
  if( ! ref( $stored_obj ) or 
      ! $stored_obj->isa( 'Bio::Root::Storable' ) ){
    $self->throw( "Token $token returned no data" );
  }
  return $stored_obj;
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

   my $SQL = qq(
SELECT object
FROM   blast_hsp
WHERE  ticket = ? );

   my $CHR_SQL = qq(
AND    chr_name = ? );

   my $RANGE_SQL = qq(
AND    chr_start <= ?
AND    chr_end   >= ? );

   my $q = $SQL;
   my @binded = ( $ticket );

   if( $chr_name ){
     $q .= $CHR_SQL;
     push @binded, $chr_name;

     if( $chr_start && $chr_end ){
       $q .= $RANGE_SQL;
       push @binded, $chr_end, $chr_start;
     }
   }
   warn( "$q: ", join( ', ',@binded ) ); 

   my $sth = $self->db->db_handle->prepare($q);
   my $rv = $sth->execute( @binded ) || $self->throw( $sth->errstr );

   my @hsps = map{ thaw( $_->[0] ) } @{$sth->fetchall_arrayref()};
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
  $self->dynamic_use( ref($hsps->[0] ) );
  my @feats = grep{ $_ } map{ $_->ens_genomic_align } @$hsps;
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

my %valid_table_types = ( HIT=>1, HSP=>1 );
sub use_date {
  my $key  = '_current_table';
  my $self = shift;
  my $type = uc( shift );
  $valid_table_types{$type} || 
    $self->throw( "Need a table type (Hit or HSP)" );

  $self->{$key} ||= {};
  if( ! $self->{$key}->{$type} ){
    my $sth = $self->db->db_handle->prepare( $SQL_SELECT_TABLE_LOG_CURRENT );
    my $rv = $sth->execute( $type ) || $self->throw( $sth->errstr );
    $rv > 0 || ( warn( "No current $type table found" ) && return );
    my $date = $sth->fetchrow_arrayref->[0];
    $date =~ s/-//g;
    $self->{$key}->{$type} = $date;
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
  my $dbh = $self->db->db_handle;

  # Get list of existing tables in database
  my $q = 'show tables';
  my $sth = $dbh->prepare( $q );
  my $rv = $sth->execute() || $self->throw( $sth->errstr );
  my $res = $sth->fetchall_arrayref;
  $sth->finish;
  my %existing_tables = map{ $_->[0], 1 } @$res;

  # Create blast_ticket, blast_result and/or blast_table_log if missing
  if( ! $existing_tables{blast_ticket} ){
    warn( "Creating blast_ticket table" );
    my $q = $SQL_CREATE_TICKET;
    my $sth = $dbh->prepare( $q );
    #my $rv = $sth->execute() || $self->throw( $sth->errstr );
  }
  if( ! $existing_tables{blast_result} ){
    warn( "Creating blast_result table" );
    my $q = $SQL_CREATE_RESULT;
    my $sth = $dbh->prepare( $q );
    #my $rv = $sth->execute() || $self->throw( $sth->errstr );    
  }
  if( ! $existing_tables{blast_table_log} ){
    warn( "Creating blast_result table" );
    my $q = $SQL_CREATE_TABLE_LOG;
    my $sth = $dbh->prepare( $q );
    #my $rv = $sth->execute() || $self->throw( $sth->errstr );    
  }

  # Get date
  my( $day, $month, $year ) = (localtime)[3,4,5];
  my $use_date = sprintf( "%04d%02d%02d", $year+1900, $month+1, $day );

  # Create today's blast_hit and blast_hsp tables if missing
  if( ! $existing_tables{'blast_hit'.$use_date} ){
    warn( "Creating todays blast_hit$use_date table" );

    # Create new table
    my $q = sprintf($SQL_CREATE_DAILY_HIT, $use_date);
    my $sth = $dbh->prepare( $q );
    my $rv = $sth->execute() || $self->throw( $sth->errstr );         
    
    # Flip current table
    my $last_use_date = $self->use_date( "HIT" ) || '';
    my $sth2 = $dbh->prepare( $SQL_TABLE_LOG_INSERT );
    my $sth3 = $dbh->prepare( $SQL_TABLE_LOG_UPDATE );
    $sth2->execute( 'CURRENT','HIT',$use_date ) 
      || die( $self->throw( $sth2->errstr ) );
    $sth3->execute( 'FILLED','HIT',$last_use_date) 
      || die( $self->throw( $sth3->errstr ) );
  }

  if( ! $existing_tables{'blast_hsp'.$use_date} ){
    warn( "Creating todays blast_hsp$use_date table" );

     # Create new table
    my $q = sprintf($SQL_CREATE_DAILY_HSP, $use_date );
    my $sth = $dbh->prepare( $q );
    my $rv = $sth->execute() || $self->throw( $sth->errstr );         
    
    my $last_use_date = $self->use_date( "HSP" ) || '';    
    my $sth2 = $dbh->prepare( $SQL_TABLE_LOG_INSERT );
    my $sth3 = $dbh->prepare(  $SQL_TABLE_LOG_UPDATE );
    $sth2->execute( 'CURRENT','HSP',$use_date ) 
      || die( $self->throw( $sth2->errstr ) );
    $sth3->execute( 'FILLED','HSP',$last_use_date) 
      || die( $self->throw( $sth3->errstr ) );
  }


#  my $days = shift || $self->throw( "Missing arg: number of days" );
#  $days =~ /\D/    && $self->throw( "Bad arg: number of days $days not int" );

#  my $q = qq/
#SELECT ticket 
#FROM   blast_ticket
#WHERE  update_time < SUBDATE( NOW(), INTERVAL $days DAY ) /;

#  my $q_del_tmpl = qq/
#DELETE
#FROM   blast_%s
#WHERE  ticket like "%s" /;

#  my $sth = $self->db->db_handle->prepare($q);
#  my $rv = $sth->execute() || $self->throw( $sth->errstr );
#  my $res = $sth->fetchall_arrayref;
#  $sth->finish;
  
#  my @types = ( 'hsp','hit','result','ticket' );
#  my %num_deleted = map{ $_=>0 } @types;

#  foreach my $row( @$res ){
#    my $ticket = $row->[0];

#    foreach my $type( @types ){
#      my $q_del = sprintf( $q_del_tmpl, $type, $ticket );
#      my $sth = $self->db->db_handle->prepare($q_del);
#      my $rv = $sth->execute() || $self->throw( $sth->errstr );
#      $num_deleted{$type} += $rv;
#    }
#  }
#  map{ warn("Purging $days days: Deleted $num_deleted{$_} rows of type $_\n") }
#  keys %num_deleted;
#  return 1;
}

#----------------------------------------------------------------------
1;

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
INSERT INTO blast_result ( ticket, object )
VALUES                   ( ? , ? )";

our $SQL_RESULT_UPDATE = "
UPDATE  blast_result
SET     object = '?'
WHERE   result_id = ?";

our $SQL_RESULT_RETRIEVE = "
SELECT object
FROM   blast_result
WHERE  result_id = ? ";

#--- HITS ---

our $SQL_HIT_STORE = "
INSERT INTO blast_hit ( ticket, object )
VALUES                ( ? , ? )";

our $SQL_HIT_UPDATE = "
UPDATE  blast_hit
SET     object = '?'
WHERE   hit_id = ?";

our $SQL_HIT_RETRIEVE = "
SELECT object
FROM   blast_hit
WHERE  hit_id = ? ";

#--- HSPS ---

our $SQL_HSP_STORE = "
INSERT INTO blast_hsp ( ticket, object, chr_name, chr_start, chr_end )
VALUES                ( ? , ? , ? , ? , ? )";

our $SQL_HSP_UPDATE = "
UPDATE  blast_hsp
SET     object    = ?,
        chr_name  = ?,
        chr_start = ?,
        chr_end   = ?
WHERE   hsp_id    = ?";

our $SQL_HSP_RETRIEVE = "
SELECT object
FROM   blast_hsp
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

  if( $rv eq '0E0' ){ # Insert (do first to minimise risk of race)
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

  my $store_obj = $res->_prepare_storable;
  my $frozen = freeze( $store_obj );
  my $ticket = $res->group_ticket;
  my $token  = $res->token;

  my( $rv );
  if( $token ){
    my $sth = $dbh->prepare( $SQL_RESULT_RETRIEVE );
    $rv = $sth->execute( $token ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }
  if( $rv ){ # Update
    my $sth = $dbh->prepare( $SQL_RESULT_UPDATE );
    $sth->execute( $store_obj ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  else{ # Insert
    my $sth = $dbh->prepare( $SQL_RESULT_STORE );
    $sth->execute( $ticket, $frozen ) || $self->throw( $sth->errstr );
    $res->token( $dbh->{mysql_insertid} );
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
  my $self   = shift;
  my $token  = shift || $self->throw( "Need a Result token" );

  my $dbh  = $self->db->db_handle;
  my $sth = $dbh->prepare( $SQL_RESULT_RETRIEVE );
  my $rv  = $sth->execute( $token ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $token not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
  $sth->finish;

  my $stored_obj = thaw( $frozen );
  if( ! ref( $stored_obj ) or 
      ! $stored_obj->isa( 'Bio::Root::Storable' ) ){
    $self->throw( "Token $token returned no data" );
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
  my $token  = $hit->token;

  my( $rv );
  if( $token ){
    my $sth = $dbh->prepare( $SQL_HIT_RETRIEVE );
    $rv = $sth->execute( $token ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }
  if( $rv ){ # Update
    my $sth = $dbh->prepare( $SQL_HIT_UPDATE );
    $sth->execute( $store_obj ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  else{ # Insert
    my $sth = $dbh->prepare( $SQL_HIT_STORE );
    $sth->execute( $ticket, $frozen ) || $self->throw( $sth->errstr );
    $hit->token( $dbh->{mysql_insertid} );
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

  my $dbh  = $self->db->db_handle;
  my $sth = $dbh->prepare( $SQL_HIT_RETRIEVE );
  my $rv  = $sth->execute( $token ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $token not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
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
  my $hsp  = shift || $self->throw( "Need an HSP obj" );
  my $dbh  = $self->db->db_handle;

  my $token = $hsp->token;

  my( $rv );
  if( $token ){
    my $sth = $dbh->prepare( $SQL_HSP_RETRIEVE );
    $rv = $sth->execute( $token ) ||  $self->throw( $sth->errstr );
    $sth->finish;
  }

  my $store_obj = $hsp->_prepare_storable;
  my $frozen = freeze( $store_obj );
  my $ticket = $hsp->group_ticket;

  my $chr_name  = 'NULL';
  my $chr_start = 'NULL';
  my $chr_end   = 'NULL';
  if( my $genomic = $hsp->genomic_hit ){
    $chr_name  = $genomic->seq_id;
    $chr_start = $genomic->start;
    $chr_end   = $genomic->end;
  } 

  if( $rv ){ # Update
    my $sth = $dbh->prepare( $SQL_HSP_UPDATE );
    my @bound = ( $frozen, $chr_name,  $chr_start, $chr_end, $token );
    $sth->execute( @bound ) || $self->throw( $sth->errstr );
    $sth->finish;
  }
  else{ # Insert
    my $sth = $dbh->prepare( $SQL_HSP_STORE );
    my @bound = ( $ticket, $frozen, $chr_name,  $chr_start, $chr_end );
    $sth->execute( @bound ) || $self->throw( $sth->errstr );
    $hsp->token( $dbh->{mysql_insertid} );
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

  my $dbh  = $self->db->db_handle;
  my $sth = $dbh->prepare( $SQL_HSP_RETRIEVE );
  my $rv  = $sth->execute( $token ) || $self->throw( $sth->errstr );
  if( $rv < 1 ){ $self->throw( "Token $token not found" ) }
  my ( $frozen ) = $sth->fetchrow_array;
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

#----------------------------------------------------------------------

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

  my $q = qq/
SELECT ticket 
FROM   blast_ticket
WHERE  update_time < SUBDATE( NOW(), INTERVAL $days DAY ) /;

  my $q_del_tmpl = qq/
DELETE
FROM   blast_%s
WHERE  ticket like "%s" /;

  my $sth = $self->db->db_handle->prepare($q);
  my $rv = $sth->execute() || $self->throw( $sth->errstr );
  my $res = $sth->fetchall_arrayref;
  $sth->finish;
  
  my @types = ( 'hsp','hit','result','ticket' );
  my %num_deleted = map{ $_=>0 } @types;

  foreach my $row( @$res ){
    my $ticket = $row->[0];

    foreach my $type( @types ){
      my $q_del = sprintf( $q_del_tmpl, $type, $ticket );
      my $sth = $self->db->db_handle->prepare($q_del);
      my $rv = $sth->execute() || $self->throw( $sth->errstr );
      $num_deleted{$type} += $rv;
    }
  }
  map{ warn("Purging $days days: Deleted $num_deleted{$_} rows of type $_\n") }
  keys %num_deleted;
  return 1;
}

#----------------------------------------------------------------------
1;

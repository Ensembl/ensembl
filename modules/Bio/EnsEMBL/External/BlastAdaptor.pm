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
  my $ticket = $hsp->blast_ticket;

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
  my $self = shift;
  my $dbh  = $self->db->db_handle;
  my( $caller, $token ) = @_;
  
  my $hsp;
  my $class = ref( $caller ) || $caller;

  # Is this a call on a retrievable object?
  if( ref( $caller ) ){
    if( $caller->retrievable ){
      $hsp = $caller;
      $token = $hsp->token;
    }
  }
  else{ $hsp = bless( {}, $caller ) }

  if( ! $token ){ $self->throw( "Need an HSP token" ) }
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
  $stored_obj->{-retrievable} = 0;
  map{ $hsp->{$_} = $stored_obj->{$_} } keys %$stored_obj; # Copy hasheys
  return bless( $hsp, ref( $stored_obj ) ); # Maintain class of stored obj

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
WHERE  \(ticket = ? OR ticket = ?\) );

   my $CHR_SQL = qq(
AND    chr_name = ? );

   my $RANGE_SQL = qq(
AND    chr_start <= ?
AND    chr_end   >= ? );

   my $q = $SQL;
   my @binded = ( $ticket, substr($ticket,6) );

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

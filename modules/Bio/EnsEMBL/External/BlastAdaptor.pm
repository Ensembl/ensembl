# Let the code begin...
package Bio::EnsEMBL::External::BlastAdaptor;

use strict;
use DBI;
use Storable qw(freeze thaw);
use Data::Dumper qw( Dumper );

use vars qw(@ISA);

use Bio::EnsEMBL::Root;

@ISA = qw( Bio::EnsEMBL::Root );

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
  my $db     = shift || $caller->throw( "Need a db_adaptor!" );
  my $self = $caller->SUPER::new();
  $self->db_adaptor($db);
  return $self;
}

#----------------------------------------------------------------------

=head2 dbh

  Arg [1]   : 
  Function  : 
  Returntype: 
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub dbh {
  my $self = shift;
  if( my $dba = shift ){
    $self->{__dbh} =$dba;
  }
  return $self->{__dbh}->db_handle;
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
AND    chr_start >= ?
AND    chr_end   <= ? );

   my $q = $SQL;
   my @binded = ( $ticket );

   if( $chr_name ){
     $q .= $CHR_SQL;
     push @binded, $chr_name;

     if( $chr_start && $chr_end ){
       $q .= $RANGE_SQL;
       push @binded, $chr_start && $chr_end;
     }
   }
   my $sth = $self->dbh->prepare($q);
   my $rv = $sth->execute( @binded ) || $self->throw( $sth->errstr );

   my @hsps = map{ thaw( $_->[0] ) } @{$sth->fetchall_arrayref()};

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
  my @feats = grep{ $_ } map{ $_->ens_genomic_align } @$hsps;
  return [ @feats ];
}




#----------------------------------------------------------------------
1;

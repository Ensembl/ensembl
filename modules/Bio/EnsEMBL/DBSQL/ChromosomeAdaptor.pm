#
# Ensembl module for Bio::EnsEMBL::DBSQL::ChromosomeAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::ChromosomeAdaptor - DB Connectivity for Chromosome Object

=head1 SYNOPSIS

$chromosome_adaptor = $db_adaptor->get_ChromosomeAdaptor();
$chromosome = $chromosome_adaptor->fetch_by_chr_name('12');

=head1 DESCRIPTION

This is a database adaptor used to retrieve chromosome objects from a database.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

=cut

package Bio::EnsEMBL::DBSQL::ChromosomeAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Chromosome;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_dbID

  Arg [1]    : int $id 
               unique database identifier for chromosome to retrieve  
  Example    : my $chromosome = $chromosome_adaptor->fetch_by_dbID(1);
  Description: Retrieves a Chromosome object from the database using its
               unique identifier.  Note the the identifier is the dbID and
               does NOT correspond to the chromosome name.  dbID 1 does NOT
               necessarily correspond to chromosome '1'
  Returntype : Bio::EnsEMBL::Chromosome
  Exceptions : thrown if $id not defined
  Caller     : general

=cut

sub fetch_by_dbID {
  my ($self,$id) = @_;

  my $chr = (); 

  unless(defined $id) {
    $self->throw("Chromosome dbID argument required\n");
  }

  unless(defined $self->{'_chr_cache'} ) {
    $self->{'_chr_cache'} = {};
    $self->{'_chr_name_cache'} = {};
  }

  #
  # If there is not already a cached version of this chromosome pull it
  # from the database and add it to the cache.
  #
  unless($chr = $self->{'_chr_cache'}->{$id}) {
    my $sth = $self->prepare( "SELECT name, length
                               FROM chromosome
                               WHERE chromosome_id = ?" );
    $sth->execute( $id );
    
    my($name, $length);
    $sth->bind_columns(\$name,\$length); 
    $sth->fetch();
   

    if($@) {
      $self->throw("Could not create chromosome from dbID $id\n" .
		   "Exception: $@\n");
    }

    unless($name) {
      $self->throw("Could determine chromosome name from dbID $id\n");
    }

    $chr = new Bio::EnsEMBL::Chromosome( -adaptor => $self,
                                         -dbID => $id,
					 -chr_name => $name,
					 -length => $length );

    $self->{'_chr_cache'}->{$id} = $chr;
    $self->{'_chr_name_cache'}->{$name} = $chr;
    
  }

  return $chr;
}


=head2 fetch_by_chr_name

  Arg [1]    : string $chr_name
               the name of the chromosome to retrieve
  Example    : $chromosome = $chromosome_adaptor->fetch_by_chr_name('X');
  Description: Retrieves a chromosome object from the database using its name.
  Returntype : Bio::EnsEMBL::Chromosome
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_chr_name{
   my ($self,$chr_name) = @_;

  my $chr = undef; 

  unless(defined $chr_name) {
    $self->throw("Chromosome name argument required\n");
  }

  unless(defined $self->{'_chr_cache'} ) {
    $self->{'_chr_cache'} = {};
    $self->{'_chr_name_cache'} = {};
  }

  #
  # If there is not already a cached version of this chromosome pull it
  # from the database and add it to the cache.
  #
  unless($chr = $self->{'_chr_name_cache'}->{$chr_name}) {
    my $sth = $self->prepare( "SELECT chromosome_id, length
                               FROM chromosome
                               WHERE name = ?" );
    $sth->execute( $chr_name );
    
    my($dbID, $length);
    $sth->bind_columns(\$dbID,\$length); 

    if ($sth->rows > 0) {
    $sth->fetch();

    $chr = new Bio::EnsEMBL::Chromosome( -adaptor => $self,
					 -dbID => $dbID,
					 -chr_name => $chr_name,
					 -length => $length );

    $self->{'_chr_cache'}->{$dbID} = $chr;
    $self->{'_chr_name_cache'}->{$chr_name} = $chr;
  }
  }
  return $chr;
}


=head2 fetch_all

  Args       : none
  Example    : @chromosomes = $chromosome_adaptor->fetch_all(); 
  Description: Retrieves every chromosome object from the database.
  Returntype : listref of Bio::EnsEMBL::Chromosome
  Exceptions : none
  Caller     : general

=cut

sub fetch_all {
  my($self) = @_;
  
  my @chrs = (); 


    my $sth = $self->prepare( "SELECT chromosome_id, name, length
                               FROM chromosome" );
    $sth->execute();
    
    my($chromosome_id, $name, $length);
    $sth->bind_columns(\$chromosome_id,\$name,\$length); 

    while($sth->fetch()) {
   
    my $chr = new Bio::EnsEMBL::Chromosome( -adaptor => $self,
					 -chr_name => $name,
					 -dbID => $chromosome_id,
					 -length => $length );

    $self->{'_chr_cache'}->{$chromosome_id} = $chr;
    $self->{'_chr_name_cache'}->{$name} = $chr;
    push @chrs, $chr;
  }

  return \@chrs;
}



sub store{
  my ($self, $chromosome) = @_;

  my $chr_name = $chromosome->chr_name;
  if(!$chr_name){
    $self->throw("can't store a chromosome without a name");
  }
  my $length = $chromosome->length;
  if(!$length){
    $self->throw("can't store a chromosome without a length");
  }
  my $sql = "insert into chromosome(name, length) values(?, ?)";

  my $sth = $self->db->prepare($sql);
  $sth->execute($chr_name, $length);
  $chromosome->dbID($sth->{'mysql_insertid'});
  $chromosome->adaptor($self);
}

1;

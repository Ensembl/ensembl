#
# Ensembl module for Bio::EnsEMBL::Lite::ChromosomeAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Lite::ChromosomeAdaptor - DB Connectivity for Chromosome Object

=head1 SYNOPSIS

$chromosome_adaptor = $db_adaptor->get_ChromosomeAdaptor();
$chromosome = $chromosome_adaptor->fetch_by_chr_name('12');

=head1 DESCRIPTION

This is a database adaptor used to retrieve chromosome objects from a database.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

=cut

package Bio::EnsEMBL::Lite::ChromosomeAdaptor;
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

  $self->throw("Chromosome dbID argument required\n") unless defined $id;

  unless(defined $self->{'_chr_cache'} ) {
    $self->{'_chr_cache'} = {};
    $self->{'_chr_name_cache'} = {};
  }

  # If there is not already a cached version of this chromosome pull it
  #     from the database and add it to the cache.
  my $chr;
  unless( $chr = $self->{'_chr_cache'}->{$id} ) {
    my $sth;
    eval {
      $sth = $self->prepare(qq(
  SELECT c.chromosome_id, c.name, c.length, cs.type, cs.value
    FROM chromosome as c left join chromosome_stat as cs on cs.chromosome_id = c.chromosome_id
   WHERE chromosome_id = ?) );
      $sth->execute( $id );
    };
    $self->throw("Could not create chromosome from dbID $id\nException: $@\n") if $@;
    my $a = $sth->fetchrow_arrayref();
    $self->throw("Could determine chromosome name from dbID $id\n") unless $a;
    $chr = $self->_create_object_from_arrayref( $a );
    if($a->[3]) {
      $chr->stat($a->[3],$a->[4]);
      while($a = $sth->fetchrow_arrayref()) { $chr->stat($a->[3],$a->[4]); }
    }
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

  unless(defined $self->{'_chr_cache'} ) {
    $self->{'_chr_cache'} = {};
    $self->{'_chr_name_cache'} = {};
  }

  my $chr;
  unless( $chr = $self->{'_chr_name_cache'}->{$chr_name}) {
    my $sth;
    eval {
      $sth = $self->prepare(qq(
  SELECT c.chromosome_id, c.name, c.length, cs.type, cs.value
    FROM chromosome as c left join chromosome_stat as cs on cs.chromosome_id = c.chromosome_id
   WHERE name = ?) );
      $sth->execute( $chr_name );
    };
    $self->throw("Could not create chromosome from chr $chr_name\nException: $@\n") if $@;
    my $a = $sth->fetchrow_arrayref();
    $self->throw("Do not recognise chromosome $chr_name\n") unless $a;
    $chr = $self->_create_object_from_arrayref( $a );
    if($a->[3]) {
      $chr->stat($a->[3], $a->[4] );
      while($a = $sth->fetchrow_arrayref()) { $chr->stat($a->[3], $a->[4]); }
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

  my $sth = $self->prepare(qq(
  SELECT c.chromosome_id, c.name, c.length, cs.type, cs.value
    FROM chromosome as c left join chromosome_stat as cs on cs.chromosome_id = c.chromosome_id) );
     $sth->execute();
  my $chr;
  my $old_id = 0;
  while( my $a = $sth->fetchrow_arrayref() ) {
    if($a->[0]==$old_id) {
      $chr->stat($a->[3], $a->[4] );
    } else {
      push @chrs, $chr if $old_id ;
      $old_id = $a->[0]; 
      $chr = $self->_create_object_from_arrayref( $a );
      if($a->[3]) {
        $chr->stat($a->[3], $a->[4]);
      }
    }
  }
  push @chrs, $chr if $old_id;
  return \@chrs;
}

=head2 _create_object_from_arrayref

  Args       : hash ref containing a row of the chromosome table
  Example    : $self->_create_object_from_hashref
  Description: Creates object from hash reference
  Returntype : Bio::Chromosome object
  Exceptions : none
  Caller     : general

=cut

sub _create_object_from_arrayref {
  my( $self,$a ) =@_;
  my $chr = new Bio::EnsEMBL::Chromosome(
    -adaptor  => $self,   -dbID   => $a->[0],
    -chr_name => $a->[1], -length => $a->[2]
  );
  return $self->{'_chr_cache'     }->{$a->[0]} = 
         $self->{'_chr_name_cache'}->{$a->[1]} = $chr ;
}

sub store{
  my ($self, $chromosome) = @_;

  $self->throw("can't store a chromosome without a name") unless my $chr_name = $chromosome->chr_name;
  $self->throw("can't store a chromosome without a length") unless my $length = $chromosome->length;
  my $sth = $self->db->prepare("insert into chromosome set name=?, length=?");
  $sth->execute($chr_name, $length);
  $chromosome->dbID($sth->{'mysql_insertid'});
  $chromosome->adaptor($self);
}
1;

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

Bio::EnsEMBL::Lite::ChromosomeAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

$chromosome_adaptor = $db_adaptor->get_ChromosomeAdaptor();
$chromosome = $chromosome_adaptor->fetch_by_chr_name('12');

=head1 DESCRIPTION

This is a database adaptor used to retrieve chromosome objects from a database.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Lite::ChromosomeAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Chromosome;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# new inherited from BaseAdaptor

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
  unless( $self->{'_chr_cache'}->{$id} ) {
    my $sth;
    eval {
      $sth = $self->prepare( "SELECT * FROM chromosome WHERE chromosome_id = ?" );
      $sth->execute( $id );
    };
    $self->throw("Could not create chromosome from dbID $id\nException: $@\n") if $@;
    my $h = $sth->fetchrow_hashref();
    $self->throw("Could determine chromosome name from dbID $id\n") unless $h;
    $self->_create_object_from_hashref( $h );
  }

  return $self->{'_chr_cache'}->{$id};
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

  unless($self->{'_chr_name_cache'}->{$chr_name}) {
    my $sth;
    eval {
      $sth = $self->prepare( "SELECT * FROM chromosome WHERE name = ?" );
      $sth->execute( $chr_name );
    };
    $self->throw("Could not create chromosome from chr $chr_name\nException: $@\n") if $@;
    my $h = $sth->fetchrow_hashref();
    $self->throw("Do not recognise chromosome $chr_name\n") unless $h;
    $self->_create_object_from_hashref( $h );
  }
  return $self->{'_chr_name_cache'}->{$chr_name};
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
  my $sth = $self->prepare( "SELECT * from chromosome");
     $sth->execute();
  while( my $h = $sth->fetchrow_hashref() ) {
    push @chrs, $self->_create_object_from_hashref( $h );
  } 
  return \@chrs;
}

sub _create_object_from_hashref {
  my( $self,$h ) =@_;
  local $_;
  my $chr = new Bio::EnsEMBL::Chromosome(
    -adaptor  => $self,        -dbID   => $h->{'chromosome_id'},
    -chr_name => $h->{'name'}, -length => $h->{'length'},
  );
  $chr->stat( $_, $h->{$_} ) foreach
    ( grep { $_ ne 'chromosome_id' && $_ ne 'name' && $_ ne 'length' } keys %$h);
  return $self->{'_chr_cache'     }->{$h->{'chromosome_id'}} = 
         $self->{'_chr_name_cache'}->{$h->{'name'}         } = $chr ;
}

=head2 get_dbID_by_chr_name

  Arg [1]    : string $chr_name
               the name of the chromosome whose dbID is wanted.
  Example    : $dbID = $chromosome_adaptor->fetch_by_dbID('X') 
  Description: Retrieves a unique database identifier for a chromosome
               using the chromosomes name.  It is not recommended that this
               method be used externally from ChromosomeAdaptor.  It should 
               probably be private and may be made private in the future. A
               better way to obtain a dbID is:
               $dbID = $chromosome_adaptor->fetch_by_chr_name('X')->dbID();
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::ChromosomeAdaptor

=cut

sub get_dbID_by_chr_name {
  my ($self, $chr_name) = @_;

  unless (defined $self->{_chr_name_mapping}) {
    $self->{_chr_name_mapping} = {};

    #get the chromo names and ids from the database
    my $sth = $self->prepare('SELECT name, chromosome_id FROM chromosome');
    $sth->execute();
    #Construct the mapping of chromosome name to id
    while( my $a=fetchrow_arrayref() ) {
      $self->{_chr_name_mapping}->{$a->[0]} = $a->[1] 
    }
  }

  return $self->{_chr_name_mapping}->{$chr_name};
}    

=head2 get_landmark_MarkerFeatures

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use Slice::get_landmark_MarkerFeatures instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_landmark_MarkerFeatures{
  my ($self,$chr_name,$glob) = @_;
  $self->warn( qq(ChromosomeAdaptor::get_landmark_MarkerFeatures is deprecated use Slice::get_landmark_MarkerFeatures instead\n) );

  $glob = 500000  unless defined $glob;

  my $statement= qq(
    SELECT chr_start, chr_end, chr_strand, name
      FROM landmark_marker 
     WHERE chr_name = '$chr_name'
     ORDER BY chr_start
  );
   
  $statement =~ s/\s+/ /g;
   
  my ($start, $end, $strand, $name);
  my $sth = $self->prepare($statement);
     $sth->execute;
     $sth->bind_columns( undef, \$start, \$end,  \$strand, \$name);
   
  my @out;
  my $sf;
  while( $sth->fetch ) {
    next if defined $sf && $sf->end + $glob > $end && $sf->id eq $name;
    $sf = Bio::EnsEMBL::SeqFeature->new();
    $sf->start(  $start  );
    $sf->end(    $end    );
    $sf->strand( $strand );
    $sf->id(     $name   );
    push @out, $sf;
  } 
  return @out;
}


=head2 get_landmark_MarkerFeatures_old

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_landmark_MarkerFeatures_old{
   my ($self,$chr_name) = @_;
   my $glob = 1000;
   $self->throw( "Method deprecated. " );
   return ();
}

1;

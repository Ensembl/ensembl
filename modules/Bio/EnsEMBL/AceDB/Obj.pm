
#
# BioPerl module for DB::Obj
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AceDB::Obj - Object representing an instance of an EnsEMBL DB

=head1 SYNOPSIS

    $db = new Bio::EnsEMBL::AceDB::Obj( -host => 'caldy' , -port => '210000' );

    $clone = $db->get_clone('X45667');

    $contig = $db->get_Contig("dJ52N12.02793");

    $gene  = $db->get_Gene('HG45501');

    

=head1 DESCRIPTION

This object represents a database that is implemented somehow (you shouldn't
care much as long as you can get the object). From the object you can pull
out other objects by their stable identifier, such as Clone (accession number),
Exons, Genes and Transcripts. The clone gives you a DB::Clone object, from
which you can pull out associated genes and features. 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


#' Let the code begin...


package Bio::EnsEMBL::AceDB::Obj;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::AceDB::Contig;
use Bio::EnsEMBL::AceDB::Clone;
use Ace;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my ($host,$port,$debug) = $self->_rearrange([qw(HOST
					   PORT
					   DEBUG
					   )],@args);

  $host || $self->throw("Database object must have a host");
  
  if( $debug ) {
     $self->_debug($debug);
 } else {
     $self->_debug(0);
 }
  
  my $ace = my $db = Ace->connect(-host => $host,
				  -port => $port);

  if( !$ace ) {
      $self->throw("Could not connect to ace database $host,$port due to " . Ace->error() . " ");
  }

  if( $self->_debug > 3 ) {
     $self->warn("Using connection $ace");
  }
     
  $self->_db_handle($ace);

# set stuff in self from @args
  return $make; # success - we hope!
}


=head2 get_Gene

 Title   : get_Gene
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Gene{
   my ($self,@args) = @_;

   $self->throw("Not implemented yet! sorry!");

}

=head2 get_Clone

 Title   : get_Clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Clone{
   my ($self,$id) = @_;

   my $clone = new Bio::EnsEMBL::DB::Clone( -id => $id,
					-dbobj => $self );
   return $clone;
}

=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig{
   my ($self,$id) = @_;

   # FIXME: should check that this id is correct in this db.

   my $contig = new Bio::EnsEMBL::AceDB::Contig ( -dbobj => $self,
					       -id => $id );

   return $contig;
}

=head2 fetch

 Title   : fetch
 Usage   :
 Function: chains to aceperl fetch underneath...
 Example :
 Returns : 
 Args    :


=cut

sub fetch{
   my ($self,$str1,$str2) = @_;

   $self->_db_handle->fetch($str1,$str2);
}


=head2 _debug

 Title   : _debug
 Usage   : $obj->_debug($newval)
 Function: 
 Example : 
 Returns : value of _debug
 Args    : newvalue (optional)


=cut

sub _debug{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_handle'} = $value;
    }
    return $self->{'_db_handle'};

}


=head2 _db_handle

 Title   : _db_handle
 Usage   : $obj->_db_handle($newval)
 Function: 
 Example : 
 Returns : value of _db_handle
 Args    : newvalue (optional)


=cut

sub _db_handle{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_handle'} = $value;
    }
    return $self->{'_db_handle'};

}

=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub DESTROY{
   my ($obj) = @_;

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}



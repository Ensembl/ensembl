
#
# BioPerl module for Contig
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::Contig - Handle onto a database stored contig

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::Contig;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::SeqFeature::Generic;
use Bio::EnsEMBL::DB::Obj;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my ($dbobj,$id) = $self->_rearrange([qw(DBOBJ
					  ID
					  )],@args);

  $id || $self->throw("Cannot make contig db object without id");
  $dbobj || $self->throw("Cannot make contig db object without db object");
  $dbobj->isa('Bio::EnsEMBL::DB::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");

  $self->id($id);
  $self->_dbobj($dbobj);

# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 seq

 Title   : seq
 Usage   : $seq = $contig->seq();
 Function: Gets a Bio::Seq object out from the contig
 Example :
 Returns : Bio::Seq object
 Args    :


=cut

sub seq{
   my ($self) = @_;
   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select sequence from dna where contig = \"$id\"");
   my $res = $sth->execute();
   # should be a better way of doing this
   my $str = $res->fetchrow_hashref->{sequence};
   
   if( ! $str ) {
       $self->throw("No DNA sequence in contig $id");
   }

   return Bio::Seq->new ( -seq => $str, -id => $id, -type => 'Dna' );
   
}

=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures{
   my ($self,@args) = @_;
   my @array;

   my $id = $self->id();

   # make the SQL query

   my $sth = $self->_dbobj->prepare("select start,end,strand,score,analysis from feature where contig = \"$id\"");
   my $res = $sth->execute();

   while( my $rowhash = $sth->fetchrow_hashref) {
      my $out = new Bio::SeqFeature::Generic;
      $out->start($rowhash->{start});
      $out->end($rowhash->{end});
      $out->strand($rowhash->{strand});
      if( defined $rowhash->{score} ) {
	  $out->score($rowhash->{score});
      }
      $out->primary_tag($rowhash->{analysis});
      $out->source_tag('EnsEMBL');
      push(@array,$out);
  }
 
   return @array;
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'id'} = $value;
    }
    return $self->{'id'};

}

=head2 _dbobj

 Title   : _dbobj
 Usage   : $obj->_dbobj($newval)
 Function: 
 Example : 
 Returns : value of _dbobj
 Args    : newvalue (optional)


=cut

sub _dbobj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};

}



#
# BioPerl module for DB::Clone
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::Clone - Object representing one clone

=head1 SYNOPSIS

    # $db is Bio::EnsEMBL::DB::Obj 

    @contig = $db->get_Contigs();

    $clone = $db->get_Clone();

    @genes    = $clone->get_all_Genes();

=head1 DESCRIPTION

Represents information on one Clone

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::Clone;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DBSQL::Contig;
use Bio::EnsEMBL::DB::CloneI;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::CloneI);

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
  $dbobj->isa('Bio::EnsEMBL::DBSQL::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");

  $self->id($id);
  $self->_dbobj($dbobj);

# set stuff in self from @args
  return $make; # success - we hope!
}


=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes{
   my ($self,@args) = @_;
   my @out;
   my $id = $self->id();
   my %got;
   # prepare the SQL statement


   my $sth = $self->_dbobj->prepare("select p3.gene from contig as p4, transcript as p3, exon_transcript as p1, exon as p2 where p4.clone = '$id' and p2.contig = p4.id and p1.exon = p2.id and p3.id = p1.transcript");

   my $res = $sth->execute();
   while( my $rowhash = $sth->fetchrow_hashref) {

       if( $got{$rowhash->{'gene'}} != 1 ) {
          my $gene = $self->_dbobj->get_Gene($rowhash->{'gene'});
	  push(@out,$gene);
	  $got{$rowhash->{'gene'}} = 1;
       }

   }
   

   return @out;


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
   my ($self,$contigid) = @_;

   # should check this contig is in this clone?
   return $self->_dbobj->get_Contig($contigid);
}

=head2 get_all_Contigs

 Title   : get_Contigs
 Usage   : foreach $contig ( $clone->get_all_Contigs ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Contigs{
   my ($self) = @_;
   my $sth;
   my @res;
   my $name = $self->id();

   my $sql = "select id from contig where clone = \"$name\" ";
#   print STDERR "Looking at $name. [$sql]\n";
   $sth= $self->_dbobj->prepare($sql);
   my $res = $sth->execute();
   my $seen = 0;

   my $count = 0;
   my $total = 0;

   while( my $rowhash = $sth->fetchrow_hashref) {
       my $contig = new Bio::EnsEMBL::DBSQL::Contig ( -dbobj => $self->_dbobj,
						   -id => $rowhash->{'id'} );
       $contig->order($count++);
       $contig->offset($total);
       
       $total += $contig->length();
       $total += 400;

       push(@res,$contig);
       $seen = 1;
   }
   if( $seen == 0  ) {
       $self->throw("Clone $name has no contigs in the database. Should be impossible, but clearly isnt...");
   }

   return @res;   
}

=head2 htg_phase

 Title   : htg_phase
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub htg_phase{
   my ($self) = @_;

   my $self = shift;
   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select htg_phase from clone where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'htg_phase'};
}

=head2 sv

 Title   : sv
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub sv{
   my ($self) = @_;

   my $self = shift;
   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select sv from clone where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'sv'};
}

=head2 embl_id

 Title   : embl_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub embl_id{
   my ($self) = @_;

   my $self = shift;
   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select embl_id from clone where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'embl_id'};
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_clone_id'} = $value;
    }
    return $obj->{'_clone_id'};

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
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_dbobj'} = $value;
    }
    return $obj->{'_dbobj'};

}

1;

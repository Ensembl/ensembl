#
# BioPerl module for Bio::EnsEMBL::TimDB::Clone
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::TimDB::Clone - Perl wrapper over Tim's directories for Clones

=head1 SYNOPSIS

    $clone = Bio::EnsEMBL::TimDB::Clone->new();
 
    $clone->add_Contig($contig);
    

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::TimDB::Clone;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::DB::CloneI;
use Bio::EnsEMBL::TimDB::Contig;

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::CloneI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called
sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize(@args);

  # set stuff in self from @args
  my ($dbobj,$id,$cgp)=$self->_rearrange([qw(DBOBJ
					     ID
					     CGP
					     )],@args);
  $id || $self->throw("Cannot make contig db object without id");
  $dbobj || $self->throw("Cannot make contig db object without db object");
  $dbobj->isa('Bio::EnsEMBL::TimDB::Obj') || 
      $self->throw("Cannot make contig db object with a $dbobj object");
  $cgp || $self->throw("Cannot make a contig db object without location data");

  # id of clone
  $self->id($id);
  # db object
  $self->_dbobj($dbobj);

  # construct and test the directory of the clone
  my $clone_dir=$dbobj->{'_unfinished_root'}."/$cgp/data/$id";
  unless(-d $clone_dir){
      $self->throw("Cannot find directory for $id");
  }
  $self->{'_clone_dir'}=$clone_dir;

  # check for sequence file
  if(!-e "$clone_dir/$id.seq"){
      $self->throw("Error: no sequence file for clone $id");
  }

  # build list of contigs for clone
  # (methods get_all_Contigs and get_Contig look at this list of objects, 
  # rather than build it)
  # FIXME
  # THERE IS CURRENTLY NO QUICK WAY TO LOOK THIS UP!!!
  # getting a list of contigs is currently horrible.  There is a list in the
  # relevant dbm file, however this has to be stepped though.  There
  # is a list in the relevant clone.seq file (headers) except when this file is 
  # absent.  There should be a .gs file for each contig.

  # open dbm for this clone (to check contigs are listed there)
  my $contig_dbm_file=$dbobj->{'_unfinished_root'}."/$cgp/unfinished_ana.dbm";
  my %unfin_contig;
  unless(dbmopen(%unfin_contig,$contig_dbm_file,0666)){
      $self->throw("Error opening contig dbm file");
  }
  my($key,$val);
  my @contig_id;
  while(($key,$val)=each %unfin_contig){
      if($key=~/^$id/){
	  push(@contig_id,$key);
	  # check for gs file
	  if(!-e "$clone_dir/$key.gs"){
	      $self->throw("Error: no gs file for contig $key");
	  }
      }
  }
  my @res;
  foreach my $contig_id (@contig_id){
      my $contig = new Bio::EnsEMBL::TimDB::Contig ( -dbobj => $self->_dbobj,
						     -id => $contig_id );
      push(@res,$contig);
  }
  $self->{'_contig_array'}=\@res;

  # build list of genes for clone
  # (method get_all_Genes look at this list of objects, rather than build it)
  # FIXME
  # THERE IS CURRENTLY NO QUICK WAY TO LOOK THIS UP!!!
  

  return $make; # success - we hope!
}


=head2 get_all_Contigs

 Title   : get_all_Contigs
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_all_Contigs{
   my ($self) = @_;
   return @{$self->{'_contig_array'}};
}


=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_Contig {
   my ($self,$contigid) = @_;
   my $c;
   foreach my $contig (@{$self->{'_contig_array'}}){
       if( $contigid eq $contig->id()){
	   $c=$contig;
	   last;
       }
   }
   unless($c){
       $self->throw("contigid $contigid not found in this clone");
   }
   return $c;
}


=head2 add_Contig

 Title   : add_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Contig{
   my ($self,$contig) = @_;
   $self->throw("Cannot add items to TimDB");
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
   my ($self) = @_;
   my %h;

   $self->throw("Tim has not reimplemented this function");

   # read into a hash to make unique
   foreach my $gene ( $self->get_all_Contigs ) {
       $h{$gene->id()} = $gene;
   }

   return values %h;

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

sub _dbobj {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_dbobj'} = $value;
    }
    return $obj->{'_dbobj'};
}

1;

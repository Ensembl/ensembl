
#
# BioPerl module for Bio::EnsEMBL::PerlDB::Clone
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::PerlDB::Clone - Pure Perl implementation of Clone

=head1 SYNOPSIS

    $clone = Bio::EnsEMBL::PerlDB::Clone->new();
 
    $clone->add_Contig($contig);
    

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::PerlDB::Clone;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::DB::CloneI;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::CloneI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'_contig_hash'} = {};
# set stuff in self from @args
 return $make; # success - we hope!
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

   return $self->{'_contig_hash'}->{$id};
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

   return values %{$self->{'_contig_hash'}};
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


   if( ! $contig->isa('Bio::EnsEMBL::DB::ContigI') ) {
       $self->warn("$contig is not a contigI object...");
   }

   $self->{'_contig_hash'}->{$contig->id()} = $contig;
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

   # read into a hash to make unique
   foreach my $contig ( $self->get_all_Contigs ) {
       foreach my $gene ( $contig->get_all_Genes ) {
	   $h{$gene->id()} = $gene;
       }

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

sub id{
   my ($obj,$value) = @_;
   if( defined $value) {
       $obj->{'id'} = $value;
   }
   return $obj->{'id'};
   
}


1;











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

  $self->{'_contig_array'} = [];
# set stuff in self from @args
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

   push(@{$self->{'_contig_array'}},$contig);
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
   foreach my $gene ( $self->get_all_Contigs ) {
       $h{$gene->id()} = $gene;
   }

   return values %h;

}


1;

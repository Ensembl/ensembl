
#
# BioPerl module for VirtualMap
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

VirtualMap - DESCRIPTION of Object

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


package Bio::EnsEMBL::DB::VirtualMap;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DB::MapContig;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;
  $self->{'mapcontighash'}= {};

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 get_MapContig

 Title   : get_MapContig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_MapContig{
   my ($self,$id) = @_;

   $id || $self->throw("Need to give an id to get a mapcontig!\n");
   
   if (my $mapcontig = $self->{'mapcontighash'}->{$id}) {
       return $mapcontig;
   }
   else {
       $self->throw("Could not find mapcontig $id");
   }
}

=head2 add_MapContig

 Title   : add_MapContig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub create_MapContig{
   my ($self,$start,$start_in,$ori,$contig) = @_;

   $start || $self->throw("Need to give a start for the MapContig!");
   $start_in || $self->throw("Need to give a start_in for the MapContig!");
   $ori || $self->throw("Need to give an orientation for the MapContig!");

   if( ! $contig->isa("Bio::EnsEMBL::DB::RawContigI") ) {
       $self->throw("$contig is not a Bio::EnsEMBL::DB::RawContig!");
   }
   
   my $id=$contig->id;
   my $mapcontig=Bio::EnsEMBL::DB::MapContig->new( -contig =>$contig,
						   -ori => $ori,
						   -start => $start,
						   -startin => $start_in
						    );
   $self->{'mapcontighash'}->{$id}=$mapcontig;
}

=head2 get_all_MapContigs

 Title   : get_all_MapContigs
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_MapContigs{
   my ($self,$reverse) = @_;
   
   my @mapcontigs = values %{$self->{'mapcontighash'}};
   if ($reverse) {
       @mapcontigs = sort { $b->start <=> $a->start} @mapcontigs;
   }
   else {
       @mapcontigs = sort { $a->start <=> $b->start} @mapcontigs;
   }
   return (@mapcontigs);
}

=head2 get_all_RawContigs

 Title   : get_all_RawContigs
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_RawContigs{
   my ($self) = @_;
   
   my @contigs;

   foreach my $mc ($self->get_all_MapContigs){
       push @contigs,$mc->contig;
   }
   return (@contigs);
}

=head2 get_all_RawContig_ids

 Title   : get_all_Raw_Contig_ids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub RawContig_ids{
   my ($self) = @_;
   
   return keys %{$self->{'mapcontighash'}};
}



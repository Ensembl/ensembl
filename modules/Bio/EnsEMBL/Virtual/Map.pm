
#
# Ensembl module for Bio::EnsEMBL::Virtual::Map
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Virtual::Map - Map of MapContigs which define a VirtualContig

=head1 SYNOPSIS

    # Virtual::Map objects are internal to VirtualContigs


=head1 DESCRIPTION

This object is basically a hash of map contigs which make a virtual contig

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Virtual::Map;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Virtual::MapContig;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

# new() is written here 

sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;
  $self->{'_contig_map'} = {};
      
# set stuff in self from @args
  return $self;
}

=head2 create_MapContig

 Title   : create_MapContig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub create_MapContig{
   my ($self,$rawcontig,$start,$end,$start_in_rawcontig,$orientation) = @_;

   if( !defined $orientation || $start > $end ||
       !ref $rawcontig || !$rawcontig->isa('Bio::EnsEMBL::DB::RawContigI') ) {
       $self->throw("Invalid arguments passed into create_MapContig ($rawcontig)");
   }

   my $out = Bio::EnsEMBL::Virtual::MapContig->new(
				       -rawcontig => $rawcontig,
				       -start => $start,
				       -end => $end,
				       -rawcontig_start => $start_in_rawcontig,
				       -orientation => $orientation
				       );

   $self->_add_MapContig($out);

   return $out;

}

=head2 get_MapContig_by_id

 Title   : get_MapContig_by_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_MapContig_by_id {
   my ($self,$name) = @_;

   return $self->{'_contig_map'}->{$name};
}

=head2 each_MapContig

 Title   : each_MapContig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_MapContig{
   my ($self,@args) = @_;

   return values %{$self->{'_contig_map'}};
}

=head2 _add_MapContig

 Title   : _add_MapContig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _add_MapContig{
   my ($self,$mc) = @_;

   $self->{'_contig_map'}->{$mc->contig->id} = $mc;
}


1;



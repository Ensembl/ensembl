
#
# BioPerl module for Bio::EnsEMBL::PerlDB::Obj
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::PerlDB::Obj - Inmemory Perl based database object for ensembl.

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


package Bio::EnsEMBL::PerlDB::Obj;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::RootI;
use Bio::EnsEMBL::DB::ObjI;
use Bio::EnsEMBL::PerlDB::Contig;
use Bio::EnsEMBL::PerlDB::Clone;


@ISA = qw(Bio::Root::RootI Bio::EnsEMBL::DB::ObjI);

sub new {
  my($class,@args) = @_;

  my $self = bless {
      _gene_hash   => {},
      _contig_hash => {},
  }, $class;

  return $self;
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
   my ($self,$geneid) = @_;

   $self->{'_gene_hash'}->{$geneid} || $self->throw("No gene with $geneid stored in this in-memory PerlDB");
   return $self->{'_gene_hash'}->{$geneid};
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
   my ($self) = @_;

   # yikes
   $self->throw("Not implemented in the object!");
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
    my ($self,$contigid)= @_;

    $self->{'_contig_hash'}->{$contigid} || $self->throw("No contig with $contigid stored in this in-memory PerlDB");
    return $self->{'_contig_hash'}->{$contigid};
}


=head2 write_Gene

 Title   : write_Gene
 Usage   : $obj->write_Gene($gene)
 Function: writes a particular gene into the database
           
 Example :
 Returns : 
 Args    :


=cut

sub write_Gene{
   my ($self,$gene) = @_;

   if( !$gene->isa("Bio::EnsEMBL::Gene") ) {
       $self->throw("$gene is not an ensembl gene. Can't store!");
   }

   $self->{'_gene_hash'}->{$gene->id()} = $gene; 
}

=head2 write_Contig

 Title   : write_Contig
 Usage   : $obj->write_Contig($contigid,$dna)
 Function: writes a contig and its dna into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Contig {
   my ($self,$contig) = @_;

   if( !$contig->isa("Bio::EnsEMBL::DB::ContigI") ) {
       $self->throw("$contig is not an ensembl contig. Can't store!");
   }

   $self->{'_contig_hash'}->{$contig->id()} = $contig; 
}



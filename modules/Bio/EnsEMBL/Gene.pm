
#
# BioPerl module for Gene
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Gene - Object for confirmed Genes

=head1 SYNOPSIS

Confirmed genes. Basically has a set of transcripts

=head1 DESCRIPTION

Needs more description.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Gene;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::SeqFeature::Generic

use Bio::Root::Object;
use Bio::EnsEMBL::Transcript;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'_transcript_array'} = [];
# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 each_unique_Exon

 Title   : each_unique_Exon
 Usage   : foreach my $exon ( $gene->each_unique_Exon )
 Function: retrieves an array of exons associated with this
           gene, made nonredudant
 Example :
 Returns : 
 Args    :


=cut

sub each_unique_Exon{
   my ($self) = @_;
   my %h;

   foreach my $trans ( $self->each_Transcript ) {
       foreach my $exon ( $trans->each_Exon ) {
	   $h{$exon->id()} = $exon;
       }
   }

   return values %h;
   
}

=head2 unique_contig_ids

 Title   : unique_contig_ids
 Usage   : foreach $id ( $gene->unique_contig_ids ) 
 Function: returns an array of contig ids made unique linked
           to this gene
 Example :
 Returns : 
 Args    :


=cut

sub unique_contig_ids{
   my ($self) = @_;
   my %h;

   foreach my $exon ( $self->each_unique_Exon() ) {
       $h{$exon->contig_id()} = 1;
   }

   return keys %h;
}


=head2 add_Transcript

 Title   : add_Transcript
 Usage   : $gene->add_Transcript($tr)
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Transcript{
   my ($self,$trans) = @_;

   if( ! $trans->isa("Bio::EnsEMBL::Transcript") ) {
       $self->throw("$trans is not a Bio::EnsEMBL::Exon!");
   }

   # at the moment, use the SeqFeature sub hash. But in the future,
   # possibly do something better?

   push(@{$self->{'_transcript_array'}},$trans);
}

=head2 each_Transcript

 Title   : each_Transcript
 Usage   : foreach $trans ( $gene->each_Transcript)
 Function:
 Example :
 Returns : An array of Transcript objects
 Args    :


=cut

sub each_Transcript {
   my ($self) = @_;
   my @sub;
   my @ret;

   return @{$self->{'_transcript_array'}};   

}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'id'} = $value;
    }
    return $obj->{'id'};

}

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: 
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}

=head2 _dump

 Title   : _dump
 Usage   : dump data structure for debugging
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _dump{
   my ($self,$fh) = @_;

   if( ! $fh ) {
       $fh = \*STDOUT;
   }

   print $fh "Gene ", $self->id(), "\n";
   foreach my $t ( $self->each_Transcript ) {
       print $fh "  Trans ", $t->id(), " :";
       foreach my $e ( $t->each_Exon ) {
	   print $fh " ",$e->id(),",";
       }
       print "\n";
   }


}


1;


#
# Ensembl module for Bio::EnsEMBL::AlignBlockSet
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AlignBlockSet - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::AlignBlockSet;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

# new() is written here 

sub new {
    my($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;
    $self->{'_align_block'}= [];
# set stuff in self from @args
    return $self;
}

=head2 get_AlignBlocks

 Title   : get_AlignBlocks
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_AlignBlocks{
   my ($self,@args) = @_;

   return @{$self->{'_align_block'}};
}


=head2 add_AlignBlock

 Title   : add_AlignBlock
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_AlignBlock{
   my ($self,$bl) = @_;

   if( !ref $bl || !$bl->isa('Bio::EnsEMBL::AlignBlock') ) {
       $self->throw("Must add an AlignBlock [$bl]");
   }

   push(@{$self->{'_align_block'}},$bl);
}


=head2 create_from_cigar

 Title   : create_from_cigar
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub create_from_cigar{
   my ($class,$dbadaptor,$cigar) = @_;

   my $self = $class->new;
   
   if( !defined $cigar ) {
       $self->throw("Must create cigar format with dbadptor and cigar line");
   }


   #cigar: hs_est 0 479 + HSHNRNPA 694 2180 + 2299.91 M 119 N 563 M 117 N 295 M 4 I 1 M 142 N 148 M 97

   if( !($cigar =~ /cigar:\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S)\s+(\S+)\s+(.*)/) ) {
       $self->throw("cannot parse cigar lines");
   }

   my $id1     = $1;
   my $start1  = $2;
   my $end1    = $3;
   my $strand1 = $4;

   my $id2     = $5;
   my $start2  = $6;
   my $end2    = $7;
   my $strand2 = $8;
   
   my $score = $9;

   my $stateline = $10;

   my @states = split(/ /,$stateline);

   my $cursor1 = $start1;
   my $cursor2 = $start2;
   while( scalar(@states) > 0 ) {
       my $state = shift @states;
       my $size  = shift @states;

       if( $size != /^\d+$/ ) {
	   $self->throw("In parsing cigar line, hit a non sized state. Suspect something is wrong in the parsing");
       }

       if( $state eq 'M' ) {
	   # this assummes the strand here is 1.
	   if( $strand1 == '-' ) {
	       $self->throw("Assummes query sequence is always strand 1. Apologies!");
	   }

	   my $alignblock = Bio::EnsEMBL::AlignBlock->new();
	   $alignblock->align_start($cursor1);
	   $alignblock->align_end($cursor1+$size);
	   $alignblock->start($cursor2);
	   $alignblock->end($cursor2+$size);
	   if( $strand2 eq '-' ) {
	       $alignblock->strand(-1);
	   } else {
	       $alignblock->strand(+1);
	   }
	   $self->add_AlignBlock($alignblock);

	   $cursor1+= $size;
	   $cursor2+= $size; 
       } elsif ( $state eq 'D' ) {
	   $cursor1 += $size;
       } elsif ( $state eq 'N' || $state eq 'I' ) {
	   $cursor2 += $size;
       } else {
	   $self->throw("Not nice - cannot figure out state $state in cigar output");
       }
   }

   return $self;
}


1;




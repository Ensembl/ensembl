
#
# BioPerl module for DB/RawContigI.pm
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::RawContigI.pm - Abstract Interface for Sequenced Contig

=head1 SYNOPSIS

This is the abstract definition of a Contig, along with 'decorator'
functions

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::RawContigI;

use strict;
use Bio::EnsEMBL::DB:ContigI;
use vars qw(@ISA);


@ISA = qw(Bio::EnsEMBL::DB::ContigI);




=head2 created

 Title   : created
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub created{
   my ($self) = @_;
   $self->throw("Class [$self] has not implemented the created method");
}

    
=head2 seq_date

 Title   : seq_date
 Usage   : $contig->seq_date()
 Function: Gives the unix time value of the dna table created datetime field, which indicates
           the original time of the dna sequence data
 Example : $contig->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date{
    my ($self) = @_;

    $self->throw("Object did not provide the seq_date method on Contig interface!");

}

=head2 get_left_overlap

 Title   : get_left_overlap
 Usage   : $overlap_object = $contig->get_left_overlap();
 Function: Returns the overlap object of contig to the left.
           This could be undef, indicating no overlap
 Returns : A Bio::EnsEMBL::ContigOverlap object
 Args    : None

=cut

sub get_left_overlap{
   my ($self,@args) = @_;

   $self->throw("Object did not provide the get_left_overlap method on Contig interface!");
}


=head2 get_right_overlap

 Title   : get_right_overlap
 Usage   : $overlap_object = $contig->get_right_overlap();
 Function: Returns the overlap object of contig to the left.
           This could be undef, indicating no overlap
 Returns : A Bio::EnsEMBL::ContigOverlap object
 Args    : None

=cut

sub get_right_overlap{
   my ($self,@args) = @_;

   $self->throw("Object did not provide the get_left_overlap method on Contig interface!");
}



1;

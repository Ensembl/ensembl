
#
# BioPerl module for DB/ContigI.pm
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::ContigI.pm - Abstract Interface for Virtual Contigs

=head1 SYNOPSIS

This is the abstract definition of a Virtual Contig. 

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DB::VirtualContigI;


use strict;
use Bio::EnsEMBL::DB::RawContigI;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DB::RawContigI);

=head2 extend_Contig

 Title   : extend_Contig
 Usage   : $newcontig = $virtualcontig->extend_Contig(1000,-1000)
 Function: Makes a new virtualcontig which has a coordinate system
           shifted by the appropiate amount on each side. 

           If this extension is impossible extends to the maximal
           amounts. 
 Example :
 Returns : A Bio::EnsEMBL::DB::VirtualContigI
 Args    :


=cut

sub extend_Contig{
   my ($self,@args) = @_;
   $self->throw("Class [$self] has not implemented the extend Contig method");

}

=head2 extend_Contig_maximally

 Title   : extend_Contig
 Usage   : $newcontig = $virtualcontig->extend_Contig_maximally()
 Function: Makes a new virtualcontig which has a coordinate system
           of the greatest size possible in this set of overlapping
           raw contigs
 Example :
 Returns : A Bio::EnsEMBL::DB::VirtualContigI
 Args    :


=cut

sub extend_Contig_maximally {
   my ($self,@args) = @_;
   $self->throw("Class [$self] has not implemented the extend Contig maximally method");

}


1;

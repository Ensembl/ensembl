
#
# Ensembl module for Bio::EnsEMBL::DBOLD::GenomicAlignAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::GenomicAlignAdaptor - DESCRIPTION of Object

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


package Bio::EnsEMBL::DBOLD::GenomicAlignAdaptor;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBOLD::BaseAdaptor;
use Bio::EnsEMBL::GenomicAlign;

@ISA = qw(Bio::EnsEMBL::DBOLD::BaseAdaptor);

# we inheriet new


=head2 fetch_GenomicAlign_by_dbID

 Title   : fetch_GenomicAlign_by_dbID
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_GenomicAlign_by_dbID{
   my ($self,$dbid) = @_;

   return Bio::EnsEMBL::GenomicAlign->new( -align_id => $dbid, -adaptor => $self);
}




=head2 get_AlignBlockSet

 Title   : get_AlignBlockSet
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_AlignBlockSet{
   my ($self,$align_id,$row_number) = @_;

   my %contighash;

   if( !defined $row_number ) {
       $self->throw("Must get AlignBlockSet by row number");
   }

   my $sth = $self->prepare("select align_start,align_end,raw_id,raw_start,raw_end,raw_strand from genomic_align_block where align_id = $align_id and align_row = $row_number order by align_start");
   $sth->execute;

   my $alignset = Bio::EnsEMBL::AlignBlockSet->new();

   while( my $ref = $sth->fetchrow_arrayref ) {
       my($align_start,$align_end,$raw_id,$raw_start,$raw_end,$raw_strand) = @$ref;
       my $alignblock = Bio::EnsEMBL::AlignBlock->new();
       $alignblock->align_start($align_start);
       $alignblock->align_end($align_end);
       $alignblock->start($raw_start);
       $alignblock->end($raw_end);
       $alignblock->strand($raw_strand);
       
       if( ! defined $contighash{$raw_id} ) {
	   $contighash{$raw_id} = $self->db->get_Contig($raw_id);
       }

       $alignblock->raw_contig($contighash{$raw_id});
       $alignset->add_AlignBlock($alignblock);
   }

   return $alignset;
}





1;










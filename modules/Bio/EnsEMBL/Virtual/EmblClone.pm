
#
# Ensembl module for Bio::EnsEMBL::Virtual::EmblClone
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Virtual::EmblClone - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Virtual::EmblClone;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Virtual::Contig;


@ISA = qw(Bio::EnsEMBL::Virtual::Contig);

sub new {
    my($class,$clone) = @_;

    my $self = {};
    bless $self,$class;

    if( !defined $clone || !ref $clone || !$clone->isa('Bio::EnsEMBL::DB::CloneI') ) {
	$self->throw("Need a CloneI implementing object, not a [$clone] for EmblClone");
    }

    $self->_make_datastructures;

    # build vmap from clone raw contig information

    my $largest_offset =0;
    my $end;

    foreach my $contig ( $clone->get_all_Contigs ) {
	$self->_vmap->create_MapContig($contig,$contig->embl_offset,$contig->embl_offset+$contig->length-1,1,1);
	if( $contig->embl_offset > $largest_offset ) {
	    $largest_offset = $contig->embl_offset;
	    $end = $contig->embl_offset + $contig->length -1;
	}
    }

    $self->_vmap->length($end);
    $self->length($end);

    $self->dbobj($clone->adaptor->db);
    $self->htg_phase($clone->htg_phase);
    return $self;
}

=head2 htg_phase

 Title   : htg_phase
 Usage   : $obj->htg_phase($newval)
 Function: 
 Returns : value of htg_phase
 Args    : newvalue (optional)


=cut

sub htg_phase{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'htg_phase'} = $value;
    }
    return $obj->{'htg_phase'};

}

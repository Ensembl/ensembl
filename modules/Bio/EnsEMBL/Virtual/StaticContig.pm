
#
# Ensembl module for Bio::EnsEMBL::Virtual::StaticContig
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Virtual::StaticContig - DESCRIPTION of Object

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


package Bio::EnsEMBL::Virtual::StaticContig;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Virtual::Contig

@ISA = qw(Bio::EnsEMBL::Virtual::Contig);

# new() is written here 

my $static_number = 0;

sub new {
    my ($class,$global_start,$vc_start_position,$global_end,@contigs) = @_;
    
    my $self = {};
    bless $self,$class;
    $self->_make_datastructures(); # back to virtual contig



    if( scalar(@contigs) == 0 ) {
	$self->throw("Cannot build a virtual contig from no raw contigs. Probably an error in the call to get raw contigs");
    }

    # assumme sorted I guess ;)
    # waaaaaay too easy

    foreach my $rc ( @contigs ) {
	$self->_vmap->create_MapContig($rc,$rc->chr_start - $global_start+$vc_start_position,$rc->chr_end - $global_start +$vc_start_position,$rc->static_golden_start,$rc->static_golden_ori);
    }

    # needs to handle overhangs...
    my @contigs = $self->_vmap->each_MapContig;
    my $last = pop @contigs;
    $self->_vmap->length($last->end);
    $self->length($last->end);
    
    $self->id("static".$static_number++);
    
    return $self;
}

1;

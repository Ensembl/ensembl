
#
# Ensembl module for Bio::EnsEMBL::DBOLD::KaryotypeAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::KaryotypeAdaptor - DESCRIPTION of Object

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


package Bio::EnsEMBL::DBOLD::KaryotypeAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::DBOLD::BaseAdaptor;

@ISA = qw(Bio::EnsEMBL::DBOLD::BaseAdaptor);


# inheriet new from BaseAdaptor



=head2 get_band_label_by_position

 Title   : get_band_label_by_position
 Usage   : $band = $kary_adp->get_band_label_by_position('chr1',10000);
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_band_label_by_position{
    my ($self,$chr,$position) = @_;

    my $sth = $self->prepare("select band from karyotype where $position <= chr_end and $position > chr_start and chr_name = '$chr'");
    $sth->execute;
    my ($band) = $sth->fetchrow_array;

    return $band;
}









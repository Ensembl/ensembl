#
# EnsEMBL module for Bio::EnsEMBL::DBOLD::RawContigAdaptor
#
# Cared for by Imre Vastrik <vastrik@ebi.ac.uk>
#
# Copyright Imre Vastrik
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::RawContigAdaptor - MySQL database adapter class for EnsEMBL Feature Objects

=head1 SYNOPSIS



=head1 DESCRIPTION



=head1 CONTACT



=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBOLD::RawContigAdaptor;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object


use Bio::EnsEMBL::DBOLD::BaseAdaptor;
use DBI;
use Bio::EnsEMBL::DBOLD::DummyStatement;



@ISA = qw(Bio::EnsEMBL::DBOLD::BaseAdaptor);

# new() is inherited from Bio::EnsEMBL::DBOLD::BaseAdaptor

sub get_internal_id_by_id
{
    my ($self, $id) = @_;
    my $sth = $self->db->prepare
    (
         "select internal_id from contig where id = '$id'"
    );
    my $res = $sth->execute;
    if(my $rowhash = $sth->fetchrow_hashref) {
	return $rowhash->{internal_id};
    } else {
	$self->warn("Could not find contig with id $id");
    }
}


1;

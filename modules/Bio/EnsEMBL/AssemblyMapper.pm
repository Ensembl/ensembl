
#
# Ensembl module for Bio::EnsEMBL::AssemblyMapper
#
# Cared for by Arne Stabenau <stabenau@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AssemblyMapper - Handles mapping from raw contigs to assembly coordinates

=head1 SYNOPSIS

   
    $map_adaptor = $dbadaptor->get_AssmeblyMapperAdaptor();
    $mapper = $map_adaptor->fetch_by_type('UCSC');
  
    $mapper->register_region('chr1',1,100000);

    my @chr_coordlist = $mapper->map_coordinates(2,5,-1,627012,"rawcontig");
 
    my @raw_coordlist = $mapper->map_coordinates(10002,10020,1,"chr1","assembly");
    

=head1 DESCRIPTION

The AssemblyMapper is a database aware mapper which handles the raw
contig to assembly mapping. It allows mappings both from raw contigs
to assembly coordinates and from assembly coordinates back to raw
contigs. 

It extracts the mapping dynamically from the database

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::AssemblyMapper;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::Root::RootI;


@ISA = qw(Bio::Root::RootI);

sub new {
  my($class,@args) = @_;

    my $self = {};
    bless $self,$class;

  $self->throw("Waiting for implementation by Ewan - placeholder object!!!");

    return $self;
}


#
#
#


1;

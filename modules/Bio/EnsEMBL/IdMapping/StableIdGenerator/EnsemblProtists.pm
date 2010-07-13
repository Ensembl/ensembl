=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblProtists;
use strict;
use warnings;
no warnings 'uninitialized';
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use base qw(Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblBacteria);

# new generator to create new protist IDs (and also deal with existing plasmodial IDs)

sub is_valid {
  my ( $self, $stable_id ) = @_;
  my $base = $self->get_base();
  return ( $stable_id and (    $stable_id =~ /$base([A-z]{1,4})(\d{11})/
                            or $stable_id =~ /PVX.*/
                            or $stable_id =~ /PKH.*/ ) );

}

sub get_base { return 'EPr' }

1;

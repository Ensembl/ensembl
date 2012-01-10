=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::OntologyDBAdaptor

=head1 DESCRIPTION

Database adaptor for the ontology database ensembl_ontology_NN.
Mostly inheriting from Bio::EnsEMBL::DBSQL::DBAdaptor, overriding its
get_available_adaptors() method.  Not doing very much else at this
moment.

=cut

package Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;

use strict;
use warnings;

use base qw ( Bio::EnsEMBL::DBSQL::DBAdaptor );

sub get_available_adaptors {
  return {
    'GOTerm' => 'Bio::EnsEMBL::DBSQL::OntologyTermAdaptor',  #deprecated
    'SOTerm' => 'Bio::EnsEMBL::DBSQL::OntologyTermAdaptor',  #deprecated
    'OntologyTerm' => 'Bio::EnsEMBL::DBSQL::OntologyTermAdaptor' };
}

1;

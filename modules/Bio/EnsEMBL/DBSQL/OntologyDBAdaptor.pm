=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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
    'MetaContainer' => 'Bio::EnsEMBL::DBSQL::MetaContainer',
    'OntologyTerm' => 'Bio::EnsEMBL::DBSQL::OntologyTermAdaptor' };
}

1;

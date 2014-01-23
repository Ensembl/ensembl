=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.
INTO 
=cut

package Bio::EnsEMBL::DBSQL::SeqRegionSynonymAdaptor;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning stack_trace_dump);
use Bio::EnsEMBL::SeqRegionSynonym;

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


sub get_synonyms{
  my $self = shift;
  my $seq_id = shift;
  
  my @results;

  my $sth = $self->prepare("select seq_region_synonym_id, synonym, external_db_id from seq_region_synonym where seq_region_id = ?");
  $sth->bind_param(1, $seq_id, SQL_INTEGER);
  $sth->execute();
  my $dbid;
  my $alt_name;
  my $ex_db;
  $sth->bind_columns(\$dbid, \$alt_name, \$ex_db);
  while($sth->fetch()){
    push @results, Bio::EnsEMBL::SeqRegionSynonym->new(-adaptor => $self,
                                                       -synonym => $alt_name,
						       -dbID => $dbid,
                                                       -external_db_id => $ex_db,
                                                       -seq_region_id => $seq_id);
  }
  $sth->finish;

  return \@results;
}

sub store {
  my $self = shift;
  my $syn = shift;

  return if($syn->is_stored($self->db));

  if(!defined($syn->seq_region_id)){
    throw("seq_region_id is needed to store a seq_region_synoym");
  }

  my $insert_ignore = $self->insert_ignore_clause();
  my $sth = $self->prepare("${insert_ignore} INTO seq_region_synonym (seq_region_id, synonym, external_db_id) VALUES (?, ?, ?)");
  $sth->bind_param(1, $syn->seq_region_id,  SQL_INTEGER);
  $sth->bind_param(2, $syn->name         ,  SQL_VARCHAR);
  $sth->bind_param(3, $syn->external_db_id, SQL_INTEGER);
  $sth->execute;
  $syn->{'dbID'} = $self->last_insert_id('seq_region_synonym_id', undef, 'seq_region_synonym');
  $sth->finish;
}



1;

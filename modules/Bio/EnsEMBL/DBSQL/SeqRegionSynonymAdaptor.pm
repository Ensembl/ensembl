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

  my $sth = $self->prepare("INSERT IGNORE INTO seq_region_synonym (seq_region_id, synonym, external_db_id) VALUES (?, ?, ?)");
  $sth->bind_param(1, $syn->seq_region_id,  SQL_INTEGER);
  $sth->bind_param(2, $syn->name         ,  SQL_VARCHAR);
  $sth->bind_param(3, $syn->external_db_id, SQL_INTEGER);
  $sth->execute;
  $syn->{'dbID'} = $sth->{'mysql_insertid'};
  $sth->finish;
}



1;

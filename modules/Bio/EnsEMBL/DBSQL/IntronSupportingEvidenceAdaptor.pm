package Bio::EnsEMBL::DBSQL::IntronSupportingEvidenceAdaptor;

=pod

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

=head1 NAME

Bio::EnsEMBL::DBSQL::IntronSupportingEvidenceAdaptor

=head1 SYNOPSIS

  my $isea = $dba->get_IntronSupportingEvidenceAdaptor();
  my $ise = $isea->fetch_by_dbID(1);
  my $ise_array = $dfa->fetch_all();

=head1 DESCRIPTION



=head1 METHODS

=cut

use strict;
use warnings;
use base qw/Bio::EnsEMBL::DBSQL::BaseAdaptor/;

use Bio::EnsEMBL::Intron;
use Bio::EnsEMBL::IntronSupportingEvidence;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;

sub fetch_by_Intron {
  my ($self, $intron) = @_;
  assert_ref($intron, 'Bio::EnsEMBL::Intron', 'intron');
  return $self->_fetch_by_exon_ids($intron->prev_Exon()->dbID(), $intron->next_Exon()->dbID());
  return;
}

sub fetch_by_Exons {
  my ($self, $prev_exon, $next_exon) = @_;
  assert_ref($prev_exon, 'Bio::EnsEMBL::Exon', 'prev_exon');
  assert_ref($next_exon, 'Bio::EnsEMBL::Exon', 'next_exon');
  return $self->_fetch_by_exon_ids($prev_exon->dbID(), $next_exon->dbID());
}

sub _fetch_by_exon_ids {
  my ($self, $prev_exon_id, $next_exon_id) = @_;
  $self->bind_param_generic_fetch($prev_exon_id, SQL_INTEGER);
  $self->bind_param_generic_fetch($next_exon_id, SQL_INTEGER);
  my $results = $self->generic_fetch('ise.previous_exon_id =? and ise.next_exon_id =?');
  return $results->[0] if @{$results};
  return;
}

sub _tables {
  return ( [ 'intron_supporting_evidence', 'ise' ] );
}

sub _columns {
  return qw/
    ise.intron_supporting_evidence_id 
    ise.hit_name ise.score ise.score_type 
    ise.previous_exon_id ise.next_exon_id
  /;
}

sub _objs_from_sth {
  my ($self, $sth) = @_;
  my $exon_adaptor = $self->db()->get_ExonAdaptor();
  my @results;
  while ( my $row = $sth->fetchrow() ) {
    my ($id, $hit_name, $score, $score_type, $prev_exon_id, $next_exon_id) = @_;
    
    my $intron = Bio::EnsEMBL::Intron->new(
      $exon_adaptor->fetch_by_dbID($prev_exon_id),
      $exon_adaptor->fetch_by_dbID($next_exon_id)
    );
    
    my $evidence = Bio::EnsEMBL::IntronSupportEvidence->new(
      -DBID => $id,
      -ADAPTOR => $self,
      -INTRON => $intron,
      -HIT_NAME => $hit_name,
      -SCORE => $score,
      -SCORE_TYPE => $score_type
    );
  }
  return @results;
}

sub store {
  my ($self, $sf) = @_;
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportEvidence', 'intron_supporting_feature');
  
  my $sql = <<SQL;
insert into intron_supporting_evidence 
(previous_exon_id, next_exon_id, hit_name, score, score_type)
values (?,?,?,?,?) 
SQL

  my $params = [
    [$sf->intron()->prev_Exon()->dbID(), SQL_INTEGER], 
    [$sf->intron()->next_Exon()->dbID(), SQL_INTEGER],
    [$sf->hit_name(), SQL_VARCHAR],
    [$sf->score(), SQL_FLOAT],
    [$sf->score_type(), SQL_VARCHAR]
  ];
  $self->dbc()->sql_helper()->execute_update(-SQL => $sql, -PARAMS => $params, -CALLBACK => sub {
    my ( $sth, $dbh ) = @_;
    $sf->dbID($self->last_insert_id());
    return;
  });
  $sf->adaptor($self);
  
  return $sf->dbID();
}

sub update {
  my ($self, $sf) = @_;
  
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportEvidence', 'intron_supporting_feature');
  
  if (! $sf->is_stored($self->db())) {
    $self->store($sf);
    return;
  }
  
  my $sql = <<'SQL';
update intron_supporting_evidence 
previous_exon_id =?, next_exon_id=?, 
hit_name=?, score=?, score_type=?
where intron_supporting_evidence_id =?
SQL
  my $params = [
    [$sf->intron()->prev_Exon()->dbID(), SQL_INTEGER], 
    [$sf->intron()->next_Exon()->dbID(), SQL_INTEGER],
    [$sf->hit_name(), SQL_VARCHAR],
    [$sf->score(), SQL_FLOAT],
    [$sf->score_type(), SQL_VARCHAR],
    [$sf->dbID(), SQL_INTEGER]
  ];
  $self->dbc()->sql_helper()->execute_update(-SQL => $sql, -PARAMS => $params);
  return;
}

sub delete {
  my ($self, $sf) = @_;
  
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportEvidence', 'intron_supporting_feature');
  
  if (! $sf->is_stored($self->db())) {
    throw "Cannot delete the supporting evidence if it has not already been stored in this database";
  }
  
  $self->dbc()->sql_helper()->execute_update(
    -SQL => 'DELETE from intron_supporting_feature where intron_supporting_feature_id =?', 
    -PARAMS => [[$sf->dbID(), SQL_INTEGER]],
  );
  
  return;
}

1;
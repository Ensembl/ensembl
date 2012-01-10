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

Bio::EnsEMBL::DBSQL::UnconventionalTranscriptAssociationAdaptor

=head1 SYNOPSIS

  $utaa = $registry->get_adaptor( 'Human', 'Core',
    'UnconventionalTranscriptAssociation' );

  my $uta = $utaa->fetch_all_by_type('antisense');

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of
UnconventionalTranscriptAssociation objects from the database.  Most of
the implementation is in the superclass BaseFeatureAdaptor.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::UnconventionalTranscriptAssociationAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::UnconventionalTranscriptAssociation;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);




=head2 fetch_all_by_interaction_type

  Arg [1]    : String type
               the type of associations to obtain
  Example    : $utas = $utaa->fetch_all_by_type('antisense');
  Description: Obtains all unconventional transcript associations that
               have a particular interaction type.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : listREF of Bio::EnsEMBL::UnconventionalTranscriptAssociations
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub fetch_all_by_interaction_type {

  my( $self, $type) = @_;

  my $sth = $self->prepare("SELECT transcript_id, gene_id, interaction_type " .
			   "FROM unconventional_transcript_association " .
			   "WHERE interaction_type = ?");

  $sth->bind_param(1, $type, SQL_VARCHAR);
  $sth->execute();

  my $results = $self->_objs_from_sth($sth);

  $sth->finish();

  return $results;

}


=head2 fetch_all_by_gene

  Arg [1]    : String gene the gene of associations to obtain
  Arg [2]    : (optional) An interaction type; if set, only associations of this type will be returned.
  Example    : $utas = $utaa->fetch_all_by_gene($gene, 'antisense');
  Description: Obtains all unconventional transcript associations that involve
               a particular gene.
  Returntype : listREF of Bio::EnsEMBL::UnconventionalTranscriptAssociations
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub fetch_all_by_gene {

  my( $self, $gene, $type) = @_;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw('$gene must be a Bio::EnsEMBL::Gene');
  }

  my $sql = "SELECT transcript_id, gene_id, interaction_type FROM unconventional_transcript_association WHERE gene_id = ?";
  $sql .= " AND interaction_type = ?" if ($type);

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $gene->dbID(), SQL_INTEGER);
  $sth->bind_param(2, $type, SQL_VARCHAR) if ($type);

  $sth->execute();

  my $results = $self->_objs_from_sth($sth);

  $sth->finish();

  return $results;

}


=head2 fetch_all_by_transcript

  Arg [1]    : String transcript the transcript of associations to obtain
  Arg [2]    : (optional) An interaction type; if set, only associations of this type will be returned.
  Example    : $utas = $utaa->fetch_all_by_transcript($transcript, 'antisense');
  Description: Obtains all unconventional transcript associations that involve
               a particular transcript.
  Returntype : listREF of Bio::EnsEMBL::UnconventionalTranscriptAssociations
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub fetch_all_by_transcript {

  my( $self, $transcript, $type) = @_;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw('$transcript must be a Bio::EnsEMBL::Transcript');
  }

  my $sql = "SELECT transcript_id, gene_id, interaction_type FROM unconventional_transcript_association WHERE transcript_id = ?";
  $sql .= " AND interaction_type = ?" if ($type);

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $transcript->dbID(), SQL_INTEGER);
  $sth->bind_param(2, $type, SQL_VARCHAR) if ($type);

  $sth->execute();

  my $results = $self->_objs_from_sth($sth);

  $sth->finish();

  return $results;

}


=head2 store

  Arg [1]    : Bio::EnsEMBL::UnconventionalTranscriptAssociation
               the unconventional transcript association to store in the database
  Example    : $utaa_adaptor->store($uta);
  Description: stores unconventional transcript associations in the database
  Returntype : none
  Exceptions :
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub store {

  my( $self, $uta ) = @_;

  if(!ref($uta) || !$uta->isa('Bio::EnsEMBL::UnconventionalTranscriptAssociation')) {
    throw('$uta must be a Bio::EnsEMBL::UnconventionalTranscriptAssociation');
  }

  my $sth = $self->prepare(qq {INSERT into unconventional_transcript_association
			       (transcript_id, gene_id, interaction_type) VALUES (?,?,?)});

  $sth->bind_param(1, $uta->transcript()->dbID(), SQL_INTEGER);
  $sth->bind_param(2, $uta->gene()->dbID,         SQL_INTEGER);
  $sth->bind_param(3, $uta->interaction_type(),   SQL_VARCHAR);

  $sth->execute();

}



sub _objs_from_sth {

  my ($self, $sth) = @_;

  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();
  my $gene_adaptor = $self->db()->get_GeneAdaptor();

  my ($gene_id, $transcript_id, $type);
  $sth->bind_columns(\$transcript_id, \$gene_id, \$type);

  my @results;

  while($sth->fetch()) {

    my $gene = $gene_adaptor->fetch_by_dbID($gene_id);
    my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);

    my $obj = Bio::EnsEMBL::UnconventionalTranscriptAssociation->new($transcript, $gene, $type);
    push @results, $obj;

  }

  return \@results;
}

1;



#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::RegulatoryFactorAdaptor
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::RegulatoryFactorAdaptor

=head1 

$rma = $database_adaptor->get_RegulatoryFactorAdaptor();

$regulatory_factor = $rfa->fetch_by_dbID(132);
$regulatory_factor = $rfa->fetch_by_name('Factor1');
@regulatory_factors = $rfa->fetch_all_by_type("transcription_factor");

$rma->store($rm1, $rm2, $rm3);

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of RegulatoryFactor objects.

=head1 AUTHOR - Glenn Proctor

Email glenn@ebi.ac.uk

=head1 CONTACT

Post questions to the EnsEMBL developer list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::RegulatoryFactorAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RegulatoryFactor;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);

use vars qw(@ISA);
@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


=head2 fetch_by_dbID

  Arg [1]    : int $db_id
               The database identifier for the RegulatoryFactor to obtain
  Example    : $regulatory_factor = $rma->fetch_by_dbID(4);
  Description: Obtains a RegulatoryFactor object from the database via its
               primary key. 
  Returntype : Bio::EnsEMBL::RegulatoryFactor
  Exceptions : none
  Caller     : general, Bio::EnsEMBL::RegulatoryFactorAdaptor
  Status     : At Risk
             : under development

=cut

sub fetch_by_dbID {
    my( $self, $db_id ) = @_;

    my ($rc) = @{$self->_generic_fetch("regulatory_factor_id = $db_id")};

    return $rc;
}



=head2 fetch_by_name

  Arg [1]    : string $name
               the name of the regulatory factor to obtain
  Example    : $rc = $rma->fetch_by_name('Factor');
  Description: Obtains a regulatory factor from the database via its name
  Returntype : Bio::EnsEMBL::RegulatoryFactor
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub fetch_by_name {
    my( $self, $name ) = @_;

    my ($rc) = @{$self->_generic_fetch("name = '$name'")};

    return $rc;
}


=head2 fetch_all_by_type

  Arg [1]    : string $type
               the type of regulatory factor to obtain
  Example    : $rm = $rma->fetch_all_by_type('promoter');
  Description: Obtains all regulatory factors of a particular type
  Returntype : listREF of Bio::EnsEMBL::RegulatoryFactors
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub fetch_all_by_type {
    my( $self, $type) = @_;

    return $self->_generic_fetch("type = '$type'");
}



=head2 _generic_fetch

  Arg [1]    : string $where_clause
  Example    : none
  Description: PRIVATE used to create RegulatoryFactor features from an
               SQL constraint
  Returntype : listref of Bio::EnsEMBL::RegulatoryFactor objects
  Exceptions : none
  Caller     : internal
  Status     : At Risk
             : under development

=cut

sub _generic_fetch {

    my ($self, $where_clause) = @_;

    my ($regulatory_factor_id, $name, $type);

    my $sth = $self->prepare(qq{
        SELECT regulatory_factor_id, name, type
        FROM regulatory_factor
        WHERE }. $where_clause);

    $sth->execute;
    $sth->bind_columns(\$regulatory_factor_id, \$name, \$type);

    my @factors;
    while ($sth->fetch) {
        push @factors, Bio::EnsEMBL::RegulatoryFactor->new(-DBID => $regulatory_factor_id,
							   -NAME => $name,
							   -TYPE => $type,
							   -ADAPTOR => $self);
      }
    return \@factors;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::RegulatoryFactors @factors
  Example    : $rma->store(@factors);
  Description: stores a list of RegulatoryFactor objects in the database
  Returntype : none
  Exceptions : none
  Caller     : ?
  Status     : At Risk
             : under development

=cut

sub store {
  my( $self, @factors ) = @_;

  my $sth = $self->prepare("INSERT into regulatory_factor (name, type) VALUES (?,?)");

  foreach my $rm (@factors) {

    my $name = $rm->name or throw("name not set");
    my $type = $rm->type or throw("type not set");

    $sth->execute($name, $type);

    my $db_id = $sth->{'mysql_insertid'}
    or throw("Didn't get an insertid from the INSERT statement");

    $rm->dbID($db_id);
    $rm->adaptor($self);
  }
}


=head2 fetch_factors_coded_for_by_gene

  Arg [1]    : Bio::EnsEMBL::Gene 
  Example    : $rfa->fetch_factors_coded_for_by_gene($gene)
  Description: Fetches any regulatory_factors that are coded for by a particular gene
  Returntype : Listref of Bio::Ensembl::RegulatoryFactor
  Exceptions : 
  Caller     : ?
  Status     : At Risk
             : under development

=cut

sub fetch_factors_coded_for_by_gene {

  my ($self, $gene) = @_;

  if (!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw('Expected Bio::Ensembl::Gene argument not [' . ref($gene) .'].');
  }

  my @factors;

  my ($factor_id);

  my $sth = $self->db()->dbc()->prepare("SELECT regulatory_factor_id
			                 FROM regulatory_factor_coding
			                 WHERE gene_id=?");

  $sth->execute($gene->dbID());
  $sth->bind_columns(\$factor_id);

  while ($sth->fetch) {
    push @factors, $self->fetch_by_dbID($factor_id);
  }

  return \@factors;

}


=head2 fetch_factors_coded_for_by_transcript

  Arg [1]    : Bio::EnsEMBL::Transcript 
  Example    : $rfa->fetch_factors_coded_for_by_transcript($transcript)
  Description: Fetches any regulatory_factors that are coded for by a particular transcript
  Returntype : Listref of Bio::Ensembl::RegulatoryFactor
  Exceptions : 
  Caller     : ?
  Status     : At Risk
             : under development

=cut

sub fetch_factors_coded_for_by_transcript {

  my ($self, $transcript) = @_;

  if (!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw('Expected Bio::Ensembl::Transcript argument not [' . ref($transcript) .'].');
  }

  my @factors;

  my ($factor_id);

  my $sth = $self->db()->dbc()->prepare("SELECT regulatory_factor_id
			                 FROM regulatory_factor_coding
			                 WHERE transcript_id=?");

  $sth->execute($transcript->dbID());
  $sth->bind_columns(\$factor_id);

  while ($sth->fetch) {
    push @factors, $self->fetch_by_dbID($factor_id);
  }

  return \@factors;

}

1;

__END__


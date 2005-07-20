#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::RegulatoryMotifAdaptor
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::RegulatoryFactorAdaptor

=head1 

$rma = $database_adaptor->get_RegulatoryFactorAdaptor();

$regulatory_factor = $rma->fetch_by_dbID(132);
$regulatory_factor = $rma->fetch_by_name('Factor1');
@regulatory_factors = $rma->fetch_all_by_type("promoter");

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
							 -TYPE => $type);
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

1;

__END__


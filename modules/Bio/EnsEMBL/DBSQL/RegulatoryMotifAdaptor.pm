#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::RegulatoryMotifAdaptor
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::RegulatoryMotifAdaptor

=head1 

$rma = $database_adaptor->get_RegulatoryMotifAdaptor();

$regulatory_motif = $rma->fetch_by_dbID(132);
$regulatory_motif = $rma->fetch_by_name('Motif1');
@regulatory_motifs = $rma->fetch_all_by_type("promoter");

$rma->store($rm1, $rm2, $rm3);

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of RegulatoryMotif objects.

=head1 AUTHOR - Glenn Proctor

Email glenn@ebi.ac.uk

=head1 CONTACT

Post questions to the EnsEMBL developer list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::RegulatoryMotifAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RegulatoryMotif;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);

use vars qw(@ISA);
@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


=head2 fetch_by_dbID

  Arg [1]    : int $db_id
               The database identifier for the RegulatoryMotif to obtain
  Example    : $regulatory_motif = $rma->fetch_by_dbID(4);
  Description: Obtains a RegulatoryMotif object from the database via its
               primary key. 
  Returntype : Bio::EnsEMBL::RegulatoryMotif
  Exceptions : none
  Caller     : general, Bio::EnsEMBL::RegulatoryMotifAdaptor

=cut

sub fetch_by_dbID {
    my( $self, $db_id ) = @_;

    my ($rc) = @{$self->_generic_fetch("regulatory_motif_id = $db_id")};

    return $rc;
}



=head2 fetch_by_name

  Arg [1]    : string $name
               the name of the regulatory motif to obtain
  Example    : $rc = $rma->fetch_by_name('Motif');
  Description: Obtains a regulatory motif from the database via its name
  Returntype : Bio::EnsEMBL::RegulatoryMotif
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
               the type of regulatory motif to obtain
  Example    : $rm = $rma->fetch_all_by_type('promoter');
  Description: Obtains all regulatory motifs of a particular type
  Returntype : listREF of Bio::EnsEMBL::RegulatoryMotifs
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
  Description: PRIVATE used to create RegulatoryMotif features from an
               SQL constraint
  Returntype : listref of Bio::EnsEMBL::RegulatoryMotif objects
  Exceptions : none
  Caller     : internal

=cut

sub _generic_fetch {

    my ($self, $where_clause) = @_;

    my ($regulatory_motif_id, $name, $type);

    my $sth = $self->prepare(qq{
        SELECT regulatory_motif_id, name, type
        FROM regulatory_motif
        WHERE }. $where_clause);

    $sth->execute;
    $sth->bind_columns(\$regulatory_motif_id, \$name, \$type);

    my @motifs;
    while ($sth->fetch) {
        push @motifs, Bio::EnsEMBL::RegulatoryMotif->new(-DBID => $regulatory_motif_id,
							 -NAME => $name,
							 -TYPE => $type);
      }
    return \@motifs;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::RegulatoryMotifs @motifs
  Example    : $rma->store(@motifs);
  Description: stores a list of RegulatoryMotif objects in the database
  Returntype : none
  Exceptions : none
  Caller     : ?

=cut

sub store {
  my( $self, @motifs ) = @_;

  my $sth = $self->prepare("INSERT into regulatory_motif (name, type) VALUES (?,?)");

  foreach my $rm (@motifs) {

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


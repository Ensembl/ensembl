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

Bio::EnsEMBL::DBSQL::SOTermAdaptor

=head1 DESCRIPTION

A specialization of Bio::EnsEMBL::DBSQL::OntologyTermAdaptor,
specifically for Sequence ontology (SO) terms.  See the
documentation of Bio::EnsEMBL::DBSQL::OntologyTermAdaptor for
further information.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::SOTermAdaptor;

use strict;
use warnings;

use base qw( Bio::EnsEMBL::DBSQL::OntologyTermAdaptor );

=head2 new

  Arg [1]       : Bio::EnsEMBL::DBSQL::DBAdaptor
                  Argument required for parent class
                  Bio::EnsEMBL::DBSQL::BaseAdaptor.

  Description   : Creates an ontology term adaptor for SO terms.

  Example       :

    my $go_adaptor = Bio::EnsEMBL::DBSQL::SOTermAdaptor->new( $dba );

  Return type   : Bio::EnsEMBL::DBSQL::SOTermAdaptor

=cut

sub new {
  my ( $proto, $dba ) = @_;

  my $this = $proto->SUPER::new( $dba, 'SO' );

  return $this;
}

1;

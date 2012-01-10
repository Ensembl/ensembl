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

package Bio::EnsEMBL::IdMapping::StableIdGenerator::AedesAegypti;

# Package that implements incrementing and verification of Aedes aegypti
# stable IDs as used by the VectorBase project.

use strict;
use warnings;

use base qw(Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric);

sub increment_stable_id {

  # This method will increment a stable ID.  For Aedes aegypti, it will
  # pick out the numerical part of the stable ID (no matter what type of
  # stable ID it is) and increment it by one.  It will then replace the
  # numerical part by the incremented value and return the new stable
  # ID.  The parsing of the stable ID is very naive.

  my ( $self, $stable_id ) = @_;

  if ( !$self->is_valid($stable_id) ) {
    throw("Unknown or missing stable ID: $stable_id.");
  }

  $stable_id =~ /^(\D*)(\d+)(\D*)/;

  my $number_as_string = "$2";
  my $number           = $2 + 1;
  $stable_id = sprintf(
    "%s" . sprintf( "%%0%dd", length($number_as_string) ) . "%s",
    $1, $number, $3 );

  return $stable_id;
}

sub is_valid {

  # A stable ID is a valid Aedes aegypti stable ID if it begins with the
  # character string "AAEL".

  my ( $self, $stable_id ) = @_;

  if ( !( defined($stable_id) && $stable_id =~ /^AAEL/ ) ) { return 0 }

  return 1;
}

1;

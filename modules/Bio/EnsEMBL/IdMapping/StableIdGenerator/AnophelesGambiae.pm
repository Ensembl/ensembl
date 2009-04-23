package Bio::EnsEMBL::IdMapping::StableIdGenerator::AnophelesGambiae;

# Package that implements incrementing and verification of Anopheles gambiae
# stable IDs as used by the VectorBase project.
# Based on Aedes_Aegypti.pm
# Differs from Aedes in that Exon stable ids like Ennnnnn not AAEL.ennnnnn
# and gene/transcript/translation start AGAP not AAEL
# also needs to explicitly exclude old-style ids, as base will look at gene_archve and will consider ENSANGG to be a higher id than any AGAP, so will use that as starting pint fornew ids

use strict;
use warnings;

use base qw(Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric);

sub increment_stable_id {

  # This method will increment a stable ID.  For Anopheles, it will
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
  $stable_id = sprintf( "%s"
                          . sprintf( "%%0%dd", length($number_as_string)
                          )
                          . "%s",
                        $1, $number, $3 );

  return $stable_id;
}

sub is_valid {

  # A stable ID is a valid Anopheles stable ID if it begins with the
  # character string "AGAP" or (for exons) just "E"
  # must make the exon one explicityly E+digits, so that will exclude old-style ENSANG ids
  # otherwise these cause problems when initial_stable_id method checks archive tables

  my ( $self, $stable_id ) = @_;

  if ( !( defined($stable_id) && ( $stable_id =~ /^AGAP/ || $stable_id =~ /^E\d+$/ ) ) ) { return 0 }

  return 1;
}


1;

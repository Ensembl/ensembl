#!/usr/bin/env perl

#-----------------------------------------------------------------------
# Demo program for the Ensembl ontology database and API.
#
# This program fetches a set of ontology terms by name and displays them
# on the console.
#-----------------------------------------------------------------------

use strict;
use warnings;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db( '-host' => 'ensembldb.ensembl.org',
                                  '-user' => 'anonymous' );

my $pattern = '%splice_site%';

# Get an ontology term adaptor and a gene adaptor (for human).
my $adaptor =
  $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

# Fetch the terms by its accession.
my @terms = @{ $adaptor->fetch_all_by_name($pattern) };

foreach my $term (@terms) {
  printf( "Accession = %s\n\tName\t= '%s'\n",
          $term->accession(), $term->name() );
  foreach my $synonym ( @{ $term->synonyms() } ) {
    printf( "\tSynonym\t= '%s'\n", $synonym );
  }
  print("\n");
}

# $Id$

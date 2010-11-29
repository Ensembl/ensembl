#!/usr/bin/env perl

#-----------------------------------------------------------------------
# Demo program for the Ensembl ontology database and API.
#
# This program fetches a GO term and uses it to retrive genes.  The
# genes retrived will be ones that are cross-referenced with either the
# GO term itself or with any of its descendant terms (following the
# transitive relation types 'is_a' or 'part_of').
#-----------------------------------------------------------------------

use strict;
use warnings;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db( '-host' => 'ensembldb.ensembl.org',
                                  '-user' => 'anonymous' );

my $accession = 'GO:0050825';

# Get a GO term adaptor and a gene adaptor (for human).
my $go_adaptor =
  $registry->get_adaptor( 'Multi', 'Ontology', 'GOTerm' );

my $gene_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );

# Fetch a GO term by its accession.
my $term = $go_adaptor->fetch_by_accession($accession);

# Use the GO term to get a bunch of genes cross-referenced to this GO
# term or to any of its descendant terms.
my @genes = @{ $gene_adaptor->fetch_all_by_GOTerm($term) };

printf( "Genes associated with the term '%s' (%s):\n",
        $term->accession(), $term->name() );

foreach my $gene (@genes) {
  printf( "stable ID = %s, external DB = %s, external name = %s\n",
          $gene->stable_id(), $gene->external_db(),
          $gene->external_name() );
}

# $Id$

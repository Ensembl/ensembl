#!/usr/bin/perl -w

#-----------------------------------------------------------------------
# Demo program for the Ensembl ontology database and API.
#
# This program fetches a GO term and displays its ancestors (following
# the transitive relation types 'is_a' and 'part_of').  The terms are
# indented with increasing distance from the original term.
#-----------------------------------------------------------------------

use strict;
use warnings;

use Bio::EnsEMBL::Registry;

# process_term() is a recursive subroutine that, given a GO term,
# displays the accessions of the parent terms via its 'is_a' or
# 'part_of' relations.  It then does the same for the parent terms.

sub process_term {
  my ( $chart, $term, $depth ) = @_;

  $depth ||= 0;

  my $indent = ' ' x ( 2*$depth );

  foreach
    my $relation ( sort keys( %{ $chart->{ $term->accession() } } ) )
  {
    if ( $relation eq 'term' ) { next }    # Skip the term itself.

    my @parent_list = @{ $chart->{ $term->accession() }{$relation} };

    # Display the accession of this term and of its parents.
    printf( "%s%s %s: %s\n",
      $indent, $term->accession(), $relation,
      join( ',', map { $_->accession() } @parent_list ) );

    foreach my $parent_term (@parent_list) {
      process_term( $chart, $parent_term, $depth + 1 );
    }
  }

} ## end sub process_term

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
  '-host' => 'ensembldb.ensembl.org',
  '-user' => 'anonymous'
);

my $accession = 'GO:0030326';

# Get a GO term adaptor.
my $go_adaptor =
  $registry->get_adaptor( 'Multi', 'Ontology', 'GOTerm' );

my $term = $go_adaptor->fetch_by_accession($accession);

my %chart = %{ $go_adaptor->_fetch_ancestor_chart($term) };

# Call the recursive process_term() subroutine with the term just
# fetched.
process_term( \%chart, $term );


# $Id$

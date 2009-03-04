=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Collection::Exon;

=head1 NAME

Bio::EnsEMBL::Collection::Exon - Feature collection implementation for
exon features.

=head1 DESCRIPTION

=head2 Extended feature representation

An exon is represented by the basic feature representation (see
documentation of Bio::EnsEMBL::Collection) and by the following extended
feature representation:

=over 4

=item 1.

Phase

=item 2.

End phase

=item 3.

Is-current

=item 4.

Stable ID

=item 5.

Version

=item 6.

Created date

=item 7.

Modified date

=back

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');

use base( 'Bio::EnsEMBL::Collection',
          'Bio::EnsEMBL::DBSQL::ExonAdaptor' );

sub _create_feature {
  my ( $this, $feature_type, $args ) = @_;

  my $feature = $this->SUPER::_create_feature( $feature_type, $args );

  if ( !$this->lightweight() ) {
    my ( $phase, $end_phase, $is_current, $stable_id, $version,
         $created_date, $modified_date )
      = rearrange( [ 'PHASE',      'END_PHASE',
                     'IS_CURRENT', 'STABLE_ID',
                     'VERSION',    'CREATED_DATE',
                     'MODIFIED_DATE'
                   ],
                   %{$args} );

    push( @{$feature},
          $phase, $end_phase, $is_current, $stable_id, $version,
          $created_date, $modified_date );
  }

  return $feature;
}

sub _create_feature_fast {
  my ( $this, $feature_type, $args ) = @_;

  throw(   '_create_feature_fast() '
         . 'is not implemented for '
         . 'exon collections' );
}

sub _tables {
  my ($this) = @_;

  my @tables = $this->SUPER::_tables();

  if ( $this->lightweight() ) { return ( $tables[0] ) }

  return @tables;
}

sub _columns {
  my ($this) = @_;

  my @columns = $this->SUPER::_columns();

  if ( $this->lightweight() ) {
    @columns[ Bio::EnsEMBL::Collection::BASIC_SLOTS .. $#columns ] =
      map( 1, Bio::EnsEMBL::Collection::BASIC_SLOTS .. $#columns );
  }

  return @columns;
}

sub _default_where_clause {
  my ($this) = @_;

  if ( $this->lightweight() ) { return '' }

  return $this->SUPER::_default_where_clause();
}

1;

# $Id$

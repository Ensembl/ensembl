=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
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

package Bio::EnsEMBL::Collection::Gene;

=head1 NAME

Bio::EnsEMBL::Collection::Gene - Feature collection implementation
for gene features.

=head1 DESCRIPTION

=head2 Extended feature representation

A gene is represented by the basic feature representation (see
documentation of Bio::EnsEMBL::Collection) and by the following extended
feature representation:

=over 4

=item 1.

Display label (display ID of the display Xref), if defined, otherwise
undef.

=item 2.

Biotype

=item 3.

Status

=item 4.

Is-current

=item 5.

Stable ID

=item 6.

Version

=item 7.

Created date

=item 8.

Modified date

=item 9.

Description

=item 10.

Confidence.

=item 11.

External name

=item 12.

External database name

=item 13.

External status

=item 14.

Source

=back

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');

use base( 'Bio::EnsEMBL::Collection',
          'Bio::EnsEMBL::DBSQL::GeneAdaptor' );

sub _create_feature {
  my ( $this, $feature_type, $args ) = @_;

  my $feature = $this->SUPER::_create_feature( $feature_type, $args );

  if ( !$this->lightweight() ) {
    my ( $stable_id,     $version,     $external_name,
         $type,          $external_db, $external_status,
         $display_xref,  $description, $created_date,
         $modified_date, $confidence,  $biotype,
         $source,        $status,      $is_current )
      = rearrange( [ 'STABLE_ID',     'VERSION',
                     'EXTERNAL_NAME', 'TYPE',
                     'EXTERNAL_DB',   'EXTERNAL_STATUS',
                     'DISPLAY_XREF',  'DESCRIPTION',
                     'CREATED_DATE',  'MODIFIED_DATE',
                     'CONFIDENCE',    'BIOTYPE',
                     'SOURCE',        'STATUS',
                     'IS_CURRENT'
                   ],
                   %{$args} );

    push( @{$feature},
          $display_xref->display_id() || undef, $biotype || $type,
          $status,          $is_current,
          $stable_id,       $version,
          $created_date,    $modified_date,
          $description,     $confidence,
          $external_name,   $external_db,
          $external_status, $source );
  }

  return $feature;
}

sub _create_feature_fast {
  my ( $this, $feature_type, $args ) = @_;

  throw(   '_create_feature_fast() '
         . 'is not implemented for '
         . 'gene collections' );
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

# $Id$

package Bio::EnsEMBL::Collection::Transcript;

=head1 NAME

Bio::EnsEMBL::Collection::Transcript - Feature collection implementation
for transcript features.

=head1 DESCRIPTION

=head2 Extended feature representation

A transcript is represented by the basic feature representation (see
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

=back

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');

use base( 'Bio::EnsEMBL::Collection',
          'Bio::EnsEMBL::DBSQL::TranscriptAdaptor' );

sub _create_feature {
  my ( $this, $feature_type, $args ) = @_;

  my $feature = $this->SUPER::_create_feature( $feature_type, $args );

  if ( !$this->lightweight() ) {
    my ( $stable_id,    $version,         $external_name,
         $external_db,  $external_status, $display_xref,
         $created_date, $modified_date,   $description,
         $biotype,      $confidence,      $external_db,
         $status,       $is_current )
      = rearrange( [ 'STABLE_ID',       'VERSION',
                     'EXTERNAL_NAME',   'EXTERNAL_DB',
                     'EXTERNAL_STATUS', 'DISPLAY_XREF',
                     'CREATED_DATE',    'MODIFIED_DATE',
                     'DESCRIPTION',     'BIOTYPE',
                     'CONFIDENCE',      'EXTERNAL_DB',
                     'STATUS',          'IS_CURRENT'
                   ],
                   %{$args} );

    push( @{$feature},
          $display_xref->display_id() || undef, $biotype,
          $status,        $is_current,
          $stable_id,     $version,
          $created_date,  $modified_date,
          $description,   $confidence,
          $external_name, $external_db,
          $external_status );
  }

  return $feature;
}

sub _create_feature_fast {
  my ( $this, $feature_type, $args ) = @_;

  throw(   '_create_feature_fast() '
         . 'is not implemented for '
         . 'transcript collections' );
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

# $Id$

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

Stable ID

=item 4.

Version

=item 5.

Created date

=item 6.

Modified date

=item 7.

Is-current

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
    my ( $phase, $end_phase, $stable_id, $version, $created_date,
         $modified_date, $is_current )
      = rearrange( [ 'PHASE',        'END_PHASE',
                     'STABLE_ID',    'VERSION',
                     'CREATED_DATE', 'MODIFIED_DATE',
                     'IS_CURRENT'
                   ],
                   %{$args} );

    push( @{$feature},
          $phase, $end_phase, $stable_id, $version, $created_date,
          $modified_date, $is_current );
  }

  return $feature;
}

sub _create_feature_fast {
  my ( $this, $feature_type, $args ) = @_;

  my $feature =
    $this->SUPER::_create_feature_fast( $feature_type, $args );

  return $feature;
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

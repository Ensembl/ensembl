# $Id$

package Bio::EnsEMBL::Collection::Gene;

use strict;
use warnings;

use base qw( Bio::EnsEMBL::Collection );

#-----------------------------------------------------------------------
# Specialized protected methods from (super) base class
# Bio::EnsEMBL::DBSQL::BaseAdaptor
#-----------------------------------------------------------------------

sub _left_join {
  return ( [ 'xref', 'x.xref_id = g.display_xref_id' ] );
}

#-----------------------------------------------------------------------
# Specialized protected methods from base class Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _extra_tables {
  return ( [ 'gene_stable_id', 'gsi' ], [ 'xref', 'x' ] );
}

sub _extra_columns {
  return ( 'g.biotype', 'g.status', 'gsi.stable_id',
           'x.display_label' );
}

sub _extra_where_clause {
  return q(
        g.is_current  = 1
    AND g.gene_id     = gsi.gene_id
  );
}

#-----------------------------------------------------------------------
# Implemented abstract protected methods from base class
# Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _feature_table { return [ 'gene', 'g' ] }

1;

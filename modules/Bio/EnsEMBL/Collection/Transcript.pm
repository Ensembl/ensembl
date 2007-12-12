# $Id$

package Bio::EnsEMBL::Collection::Transcript;

use strict;
use warnings;

use base qw( Bio::EnsEMBL::Collection );

#-----------------------------------------------------------------------
# Specialized protected methods from base class Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _extra_tables {
  return ( [ 'transcript_stable_id', 'tsi' ],
           [ 'gene',           'g' ],
           [ 'gene_stable_id', 'gsi' ] );
}

sub _extra_columns {
  return ( 't.biotype', 't.status', 'tsi.stable_id', 'gsi.stable_id' );
}

sub _extra_where_clause {
  return q( t.is_current = 1
    AND t.transcript_id = tsi.transcript_id
    AND t.gene_id       = g.gene_id
    AND g.gene_id       = gsi.gene_id
  );
}

#-----------------------------------------------------------------------
# Implemented abstract protected methods from base class
# Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _feature_table { return [ 'transcript', 't' ] }

1;

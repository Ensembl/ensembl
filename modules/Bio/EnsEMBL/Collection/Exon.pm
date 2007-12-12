# $Id$

package Bio::EnsEMBL::Collection::Exon;

use strict;
use warnings;

use base qw( Bio::EnsEMBL::Collection );

#-----------------------------------------------------------------------
# Specialized protected methods from base class Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _extra_tables {
  return ( [ 'exon_stable_id',       'esi' ],
           [ 'exon_transcript',      'et' ],
           [ 'transcript',           't' ],
           [ 'transcript_stable_id', 'tsi' ],
           [ 'gene',                 'g' ],
           [ 'gene_stable_id',       'gsi' ] );
}

sub _extra_columns {
  return ( 'esi.stable_id', 'tsi.stable_id', 'gsi.stable_id' );
}

sub _extra_where_clause {
  return q( e.is_current = 1
    AND e.exon_id           = esi.exon_id
    AND e.exon_id           = et.exon_id
    AND et.transcript_id    = t.transcript_id
    AND t.transcript_id     = tsi.transcript_id
    AND t.gene_id           = g.gene_id
    AND g.gene_id           = gsi.gene_id
  );
}

#-----------------------------------------------------------------------
# Implemented abstract protected methods from base class
# Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _feature_table { return [ 'exon', 'e' ] }

1;

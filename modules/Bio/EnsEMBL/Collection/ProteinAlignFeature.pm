# $Id$

package Bio::EnsEMBL::Collection::ProteinAlignFeature;

use strict;
use warnings;

use base qw( Bio::EnsEMBL::Collection );

#-----------------------------------------------------------------------
# Specialized protected methods from super base class
# Bio::EnsEMBL::DBSQL::BaseAdaptor
#-----------------------------------------------------------------------

# sub _straight_join { return 1 }

#-----------------------------------------------------------------------
# Specialized protected methods from base class Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

# sub _extra_tables { }

sub _extra_columns {
  return ( 'paf.analysis_id',    'paf.hit_start',
           'paf.hit_end',        'paf.hit_name',
           'paf.cigar_line',     'paf.evalue',
           'paf.perc_ident',     'paf.score',
           'paf.external_db_id', 'paf.hcoverage' );
}

# sub _extra_where_clause { }

#-----------------------------------------------------------------------
# Implemented abstract protected methods from base class
# Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _feature_table { return [ 'protein_align_feature', 'paf' ] }

1;

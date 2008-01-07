# $Id$

package Bio::EnsEMBL::Collection::DnaAlignFeature;

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
  return ( 'daf.analysis_id', 'daf.hit_start',
           'daf.hit_end',     'daf.hit_strand',
           'daf.hit_name',    'daf.cigar_line',
           'daf.evalue',      'daf.perc_ident',
           'daf.score',       'daf.external_db_id',
           'daf.hcoverage' );
}

# sub _extra_where_clause { }

#-----------------------------------------------------------------------
# Implemented abstract protected methods from base class
# Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _feature_table { return [ 'dna_align_feature', 'daf' ] }

1;

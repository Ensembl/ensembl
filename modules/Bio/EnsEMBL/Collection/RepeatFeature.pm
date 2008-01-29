# $Id$

package Bio::EnsEMBL::Collection::RepeatFeature;

use strict;
use warnings;

use base qw( Bio::EnsEMBL::Collection );

#-----------------------------------------------------------------------
# Specialized protected methods from super base class
# Bio::EnsEMBL::DBSQL::BaseAdaptor
#-----------------------------------------------------------------------

sub _straight_join {
  # Makes no difference? The RepeatFeature adaptor has it...
  return 1;
}

#-----------------------------------------------------------------------
# Specialized protected methods from base class Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _extra_tables {
  return ( [ 'repeat_consensus', 'rc' ] );
}

sub _extra_columns {
  return ( 'rf.analysis_id', 'rf.repeat_start',
           'rf.repeat_end',  'rf.score',
           'rc.repeat_name', 'rc.repeat_class',
           'rc.repeat_type', 'rc.repeat_consensus' );
}

sub _extra_where_clause {
  return 'rc.repeat_consensus_id = rf.repeat_consensus_id';
}

# sub _has_analysis { return 1 }

#-----------------------------------------------------------------------
# Implemented abstract protected methods from base class
# Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _feature_table { return [ 'repeat_feature', 'rf' ] }

1;

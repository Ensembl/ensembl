# $Id$

package Bio::EnsEMBL::Collection::RepeatFeature;

use strict;
use warnings;

use base qw( Bio::EnsEMBL::Collection );

#-----------------------------------------------------------------------
# Specialized protected methods from base class Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _extra_tables {
  return ( [ 'repeat_consensus', 'rc' ], [ 'analysis', 'a' ] );
}

sub _extra_columns {
  return ( 'rf.repeat_start', 'rf.repeat_end',
           'rf.score',        'rc.repeat_name',
           'a.logic_name' );
}

sub _extra_where_clause {
  return q( rc.repeat_consensus_id = rf.repeat_consensus_id
    AND a.analysis_id = rf.analysis_id
  );
}

#-----------------------------------------------------------------------
# Implemented abstract protected methods from base class
# Bio::EnsEMBL::Collection
#-----------------------------------------------------------------------

sub _feature_table { return [ 'repeat_feature', 'rf' ] }

1;

use strict;
use Bio::EnsEMBL::Utils::Cache;


package Bio::EnsEMBL::Utils::SeqRegionCache;

my $SEQ_REGION_CACHE_SIZE = 4000;

my %sr_id_cache;
my %sr_name_cache;

tie(%sr_name_cache, 'Bio::EnsEMBL::Utils::Cache', $SEQ_REGION_CACHE_SIZE);
tie(%sr_name_cache, 'Bio::EnsEMBL::Utils::Cache', $SEQ_REGION_CACHE_SIZE);

1;


#
# the items to cache should be listrefs to
# [ sr_id, sr_name, cs_id, sr_length ]
#
# The name cache key is "sr_name:cs_id"
# The id cache is keyed on "sr_id"
#


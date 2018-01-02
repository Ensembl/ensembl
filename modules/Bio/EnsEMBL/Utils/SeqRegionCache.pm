=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Utils::SeqRegionCache - A shared LRU cache of information about
seq_regions

=head1 SYNOPSIS

  use Bio::EnsEMBL::DBSQL::DBAdaptor;

  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  $seq_region_cache = $db->get_SeqRegionCache();

  $key = "$seq_region_name:$coord_system_id";

  $array = $seq_region_cache->{$key};

  if ($array) {
    $name   = $array->[1];
    $length = $array->[3];
  } else {
    # cache miss, get the info from the database
    # ...

    # cache the retrieved information
    $seq_region_cache->{$key} = [
      $seq_region_id,   $seq_region_name,
      $coord_system_id, $seq_region_length
    ];
  }

=head1 DESCRIPTION

This module is simply a convenient place to put a cache of sequence
region information which is shared by several adaptors for a given
database.

=head1 METHODS

=cut

use strict;
use Bio::EnsEMBL::Utils::Cache;

package Bio::EnsEMBL::Utils::SeqRegionCache;

our $SEQ_REGION_CACHE_SIZE = 40000;



sub new {
  my $class = shift;

  my %id_cache;
  my %name_cache;

  #
  # the items to cache should be listrefs to
  # [ sr_id, sr_name, cs_id, sr_length ]
  #
  # The name cache key is "sr_name:cs_id"
  # The id cache is keyed on "sr_id"
  #

  tie(%name_cache, 'Bio::EnsEMBL::Utils::Cache', $SEQ_REGION_CACHE_SIZE);
  tie(%id_cache, 'Bio::EnsEMBL::Utils::Cache', $SEQ_REGION_CACHE_SIZE);

  return bless {'name_cache' => \%name_cache, 
                'id_cache' => \%id_cache}, $class;
}


1;





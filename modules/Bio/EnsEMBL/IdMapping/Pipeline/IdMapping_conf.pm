=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::IdMapping::Pipeline::IdMapping_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
  my ($self) = @_;

  return {
    %{$self->SUPER::default_options()},
    'hive_use_param_stack'  => 1,
    'work_dir'              => $self->o('ENV', 'BASE_DIR'),
    'bin_dir'               => $self->o('work_dir')."/ensembl/misc-scripts/id_mapping",
    'config'                => $self->o('bin_dir')."/default.conf",
    'mode'                  => 'normal',
    'logauto'               => 1,
    'no_check'              => 0,
    'cache_method'          => 'build_cache_auto',
    'drop_backup_tables'    => 1,
    'backup_tables'         => 1,
    'delete_from_tables'    => 1
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'manage_id_mapping_tables',
      -module     => 'Bio::EnsEMBL::IdMapping::Pipeline::ManageIdMappingTables',
      -input_ids  => [{}],
      -parameters => {
        'config'       => $self->o('config'),
        'drop_backup_tables' => $self->o('drop_backup_tables'),
        'backup_tables' => $self->o('backup_tables'),
        'delete_from_tables' => $self->o('delete_from_tables'),
      },
      -flow_into  => ['init_check']
    },
    {
      -logic_name => 'init_check',
      -module     => 'Bio::EnsEMBL::IdMapping::Pipeline::InitCheck',
      -parameters => {
        'config'   => $self->o('config'),
        'mode'     => $self->o('mode'),
        'logauto'  => $self->o('logauto'),
        'no_check' => $self->o('no_check')
      },
      -flow_into  => ['schedule_dump_cache']
    },
    {
      -logic_name => 'schedule_dump_cache',
      -module     => 'Bio::EnsEMBL::IdMapping::Pipeline::ScheduleDumpCache',
      -parameters => {
        'cache_method' => $self->o('cache_method'),
      },
      -flow_into  => {
        '2->A' => 'dump_cache_by_seq_region',
        'A->1' => WHEN('#mode# eq "upload"' => 'upload',
                  ELSE 'mapping')
      },
    },
    {
      -logic_name => 'dump_cache_by_seq_region',
      -module     => 'Bio::EnsEMBL::IdMapping::Pipeline::DumpCacheBySeqRegion',
    },
    {
      -logic_name => 'mapping',
      -module     => 'Bio::EnsEMBL::IdMapping::Pipeline::Mapping',
      -flow_into  => ['upload']
    },
    {
      -logic_name => 'upload',
      -module     => 'Bio::EnsEMBL::IdMapping::Pipeline::Upload',
      -flow_into  => ['analyse_results']
    },
    {
      -logic_name => 'analyse_results',
      -module     => 'Bio::EnsEMBL::IdMapping::Pipeline::AnalyseResults',
    },
  ];
}

1;


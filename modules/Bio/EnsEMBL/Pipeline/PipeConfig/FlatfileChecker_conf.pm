package Bio::EnsEMBL::Pipeline::PipeConfig::FlatfileChecker_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
    my ($self) = @_;
    
    return {
      # inherit other stuff from the base class
      %{ $self->SUPER::default_options() }, 
      
      # 'base_path' => '', #where do you want your files
      # 'type' => '',
      
      ### Defaults 
      
      pipeline_name => 'flatfile_dump_check_'.$self->o('format'),
      
      pipeline_db => {
        -driver => 'sqlite',
      }
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands}, 
    ];
}

sub pipeline_analyses {
  my ($self) = @_;
  return [
    {
      -logic_name => 'FindFiles',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -meadow_type => 'LOCAL',
      -parameters => {
        inputcmd => 'find '.$self->o('base_path').q{ -type f -name '*.dat.gz'},
        column_names => ['file'],
        randomize => 1,
        input_id => '{ file => "#file#" }'
      },
      -input_ids  => [ {} ],
      -flow_into  => {
        # 1 => 'Notify',
        2 => ['CheckFlatfile'],
      },
    },
    {
      -logic_name => 'CheckFlatfile',
      -module     => 'Bio::EnsEMBL::Pipeline::Flatfile::CheckFlatfile',
      -hive_capacity => 15,
      -rc_name => 'dump',
    },
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
    %{ $self->SUPER::pipeline_wide_parameters() },
    format => $self->o('type'),
  };
}

sub resource_classes {
  my $self = shift;
  return {
    %{$self->SUPER::resource_classes()},
    dump => { 'LSF' => '-q normal -M3000000 -R"select[mem>3000] rusage[mem=3000]"'},
  }
}

1;
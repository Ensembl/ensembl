package Bio::EnsEMBL::Pipeline::PipeConfig::Core_handover_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
    my ($self) = @_;
    
    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() }, 
        
        ### OVERRIDE
        
        ### Optional overrides        
        species => [],
        
        release => software_version(),

        run_all => 0,

        bin_count => '150',

        max_run => '100',
        
        ### Defaults 
        
        pipeline_name => 'core_handover_update_'.$self->o('release'),
        
        email => $self->o('ENV', 'USER').'@sanger.ac.uk',
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands}, 
    ];
}

## See diagram for pipeline structure 
sub pipeline_analyses {
    my ($self) = @_;
    
    return [
    
      {
        -logic_name => 'ScheduleSpecies',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::ClassSpeciesFactory',
        -parameters => {
          species => $self->o('species'),
          run_all => $self->o('run_all')
        },
        -input_ids  => [ {} ],
        -max_retry_count  => 10,
        -flow_into  => {
         '3->B'  => ['PercentRepeat'],
         'B->3'  => ['PercentGC'], 
         '3->C'  => ['CodingDensity'],
         'C->3'  => ['NonCodingDensity'],
         '5->A'  => ['ConstitutiveExonsVega', 'PepStatsVega'],
         '3->A'  => ['PercentRepeat', 'CodingDensity', 'NonCodingDensity', 'PercentGC'],
         '2->A'  => ['GeneGC', 'PepStats', 'GeneCount', 'ConstitutiveExons'],
         'A->1'  => ['Notify'], 
        },
      },

      {
        -logic_name => 'ConstitutiveExons',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::ConstitutiveExons',
        -parameters => {
          dbtype => 'core',
        },
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'normal',
      },

      {
        -logic_name => 'PepStats',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::PepStats',
        -parameters => {
          tmpdir => '/tmp', binpath => '/software/pubseq/bin/emboss',
          dbtype => 'core',
        },
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'mem',
      },

      {
        -logic_name => 'ConstitutiveExonsVega',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::ConstitutiveExons',
        -parameters => {
          dbtype => 'vega',
        },
        -max_retry_count  => 5,
        -hive_capacity    => 10,
        -rc_name          => 'normal',
        -can_be_empty     => 1,
      },

      {
        -logic_name => 'PepStatsVega',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::PepStats',
        -parameters => {
          tmpdir => '/tmp', binpath => '/software/pubseq/bin/emboss',
          dbtype => 'vega',
        },
        -max_retry_count  => 5,
        -hive_capacity    => 10,
        -rc_name          => 'mem',
        -can_be_empty     => 1,
      },

      {
        -logic_name => 'GeneCount',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::GeneCount',
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'normal',
      },

      {
        -logic_name => 'NonCodingDensity',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::NonCodingDensity',
        -parameters => {
          logic_name => 'noncodingdensity', value_type => 'sum',
        },
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'normal',
        -can_be_empty     => 1,
        -flow_into => ['PseudogeneDensity'],
      },

      {
        -logic_name => 'PseudogeneDensity',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::PseudogeneDensity',
        -parameters => {
          logic_name => 'pseudogenedensity', value_type => 'sum',
        },
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'normal',
        -can_be_empty     => 1,
      },

      {
        -logic_name => 'CodingDensity',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::CodingDensity',
        -parameters => {
          logic_name => 'codingdensity', value_type => 'sum',
        },
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'normal',
        -can_be_empty     => 1,
      },

      {
        -logic_name => 'GeneGC',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::GeneGC',
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name => 'normal',
      },

      {
        -logic_name => 'PercentGC',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::PercentGC',
        -parameters => {
          table => 'repeat', logic_name => 'percentgc', value_type => 'ratio',
        },
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'normal',
        -can_be_empty     => 1,
      },

      {
        -logic_name => 'PercentRepeat',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::PercentRepeat',
        -parameters => {
          logic_name => 'percentagerepeat', value_type => 'ratio',
        },
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'mem',
        -can_be_empty     => 1,
      },

      ####### NOTIFICATION
      
      {
        -logic_name => 'Notify',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::EmailSummaryCore',
        -parameters => {
          email   => $self->o('email'),
          subject => $self->o('pipeline_name').' has finished',
        },
      }
    
    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    
    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
        release => $self->o('release'),
        bin_count => $self->o('bin_count'),
        max_run => $self->o('max_run'),
    };
}

# override the default method, to force an automatic loading of the registry in all workers
sub beekeeper_extra_cmdline_options {
    my $self = shift;
    return "-reg_conf ".$self->o("registry");
}

sub resource_classes {
    my $self = shift;
    return {
      'default' => { 'LSF' => ''},
      'normal'  => { 'LSF' => '-q normal -M 500000 -R"select[mem>500 && myens_stag1tok>800 && myens_stag2tok>800] rusage[mem=500:myens_stag1tok=10:myens_stag2tok=10:duration=10]"'},
      'mem'     => { 'LSF' => '-q normal -M 1500000 -R"select[mem>1500 && myens_stag1tok>800 && myens_stag2tok>800] rusage[mem=1500:myens_stag1tok=10:myens_stag2tok=10:duration=10]"'},
    }
}

1;

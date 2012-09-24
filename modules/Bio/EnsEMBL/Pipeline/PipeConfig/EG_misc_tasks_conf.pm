package Bio::EnsEMBL::Pipeline::PipeConfig::EG_misc_tasks_conf;

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

        division => [],

        release => software_version(),

        run_all => 0,

        bin_count => '150',

        max_run => '100',

        ### Defaults

        pipeline_name => 'misc_tasks_'.$self->o('release'),

        email => $self->o('ENV', 'USER').'@ebi.ac.uk',
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
        -module     => 'Bio::EnsEMBL::Pipeline::Production::EGSpeciesFactory',
        -parameters => {
          species => $self->o('species'),
          division => $self->o('division'),
          run_all => $self->o('run_all'),
          max_run => $self->o('max_run')

        },
        -input_ids  => [ {} ],
        -max_retry_count  => 1,
        -flow_into  => {
         '3->B' => ['PercentRepeat'],
         'B->3' => ['PercentGC'],
         '3->C' => ['CodingDensity'],
         'C->3' => ['NonCodingDensity'],
         '3->A' => ['PercentRepeat', 'CodingDensity', 'NonCodingDensity', 'PercentGC'],
         '2->A' => ['GeneGC', 'GeneCount', 'ConstitutiveExons'], # Should inclued 'PepStats'
         'A->1' => ['NotifyCore'],
         '4->D' => ['SnpDensity', 'SnpCount'],
         'D->1' => ['NotifyVariation'],
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
          tmpdir => '/tmp', binpath => '/nfs/panda/ensemblgenomes/external/EMBOSS',
          dbtype => 'core',
        },
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'mem',
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
          bin_count => $self->o('bin_count'), max_run => $self->o('max_run'),
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
          bin_count => $self->o('bin_count'), max_run => $self->o('max_run'),
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
          bin_count => $self->o('bin_count'), max_run => $self->o('max_run'),
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
          bin_count => $self->o('bin_count'), max_run => $self->o('max_run'),
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
          bin_count => $self->o('bin_count'), max_run => $self->o('max_run'),
        },
        -max_retry_count  => 3,
        -hive_capacity    => 100,
        -rc_name          => 'mem',
        -can_be_empty     => 1,
      },

      {
        -logic_name => 'SnpCount',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::SnpCount',
        -max_retry_count  => 1,
        -hive_capacity    => 10,
        -rc_name          => 'normal',
        -can_be_empty     => 1,
      },

      {
        -logic_name => 'SnpDensity',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::SnpDensity',
        -parameters => {
          table => 'gene', logic_name => 'snpdensity', value_type => 'sum',
          bin_count => $self->o('bin_count'), max_run => $self->o('max_run'),
        },
        -max_retry_count  => 1,
        -hive_capacity    => 10,
        -rc_name          => 'normal',
        -can_be_empty     => 1,
      },

      ####### NOTIFICATION

      {
        -logic_name => 'NotifyCore',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::EmailSummaryCore',
        -parameters => {
          email   => $self->o('email'),
          subject => $self->o('pipeline_name').' (core) has finished',
        },
      },

      {
        -logic_name => 'NotifyVariation',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::EmailSummaryVariation',
        -parameters => {
          email   => $self->o('email'),
          subject => $self->o('pipeline_name').' (variation) has finished',
        },
      }

    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
        release => $self->o('release'),
        species => $self->o('species'),
        species => $self->o('division'),
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
      'normal'  => { 'LSF' => '-q production-rh6 -n 4 -M 4000 -R "rusage[mem=4000]"'},
      'mem'     => { 'LSF' => '-q production-rh6 -n 4 -M 12000 -R "rusage[mem=12000]"'},
    }
}

1;

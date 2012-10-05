package Bio::EnsEMBL::Pipeline::PipeConfig::Variation_handover_conf;

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
        
        pipeline_name => 'variation_handover_update_'.$self->o('release'),
        
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
          run_all => $self->o('run_all'),
          max_run => $self->o('max_run')
        },
        -input_ids  => [ {} ],
        -flow_into  => {
          'A->1' => ['Notify'],
          '4->A' => ['SnpDensity', 'SnpCount', 'NonSense'],
        },
      },

      {
        -logic_name => 'SnpCount',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::SnpCount',
        -max_retry_count  => 2,
        -hive_capacity    => 10,
        -rc_name          => 'default',
      },

      {
        -logic_name => 'SnpDensity',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::SnpDensity',
        -parameters => {
          table => 'gene', logic_name => 'snpdensity', value_type => 'sum',
          bin_count => $self->o('bin_count'), max_run => $self->o('max_run'),
        },
        -max_retry_count  => 2,
        -hive_capacity    => 10,
        -rc_name          => 'default',
      },

      {
        -logic_name => 'NonSense',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::NonSense',
        -parameters => {
          frequency => 0.1, observation => 20,
        },
        -max_retry_count  => 2,
        -hive_capacity    => 10,
        -rc_name          => 'default',
      },

      ####### NOTIFICATION
      
      {
        -logic_name => 'Notify',
        -module     => 'Bio::EnsEMBL::Pipeline::Production::EmailSummaryVariation',
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
      'default' => { 'LSF' => '-R"select[myens_stag1tok>800 && myens_stag2tok>800] rusage[myens_stag1tok=10:myens_stag2tok=10:duration=10]"'},
      'normal'  => { 'LSF' => '-q normal -M 500000 -R"select[mem>500 && myens_stag1tok>800 && myens_stag2tok>800] rusage[mem=500:myens_stag1tok=10:myens_stag2tok=10:duration=10]"'},
      'mem'     => { 'LSF' => '-q normal -M 1000000 -R"select[mem>1000 && myens_stag1tok>800 && myens_stag2tok>800] rusage[mem=1000:myens_stag1tok=10:myens_stag2tok=10:duration=10]"'},
    }
}

1;

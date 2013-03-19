package Bio::EnsEMBL::Pipeline::PipeConfig::NcbiBlastReIndexer_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
    my ($self) = @_;
    
    my $settings = {
      # inherit other stuff from the base class
      %{ $self->SUPER::default_options() }, 
      
      # 'registry'  => '',  # registry file still needed for all species
      # 'fa_dir'    => '',  # locate where the current fa files to look for are
      # 'base_path' => '',  # where do you want your files
      
      ncbiblast_exe => 'makeblastdb',
      
      ### Defaults 
      
      pipeline_name => 'ncbi_blast_indexer',
      
    };
    $settings->{'pipeline_db'}->{'-driver'} = 'sqlite';
    return $settings;
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
        inputcmd => 'find '.$self->o('fa_dir').q{ -type f -name '#pattern#'},
        column_names => ['file']
      },
      -input_ids  => [ 
        { pattern => '*.dna*.toplevel.fa.gz' },
        { pattern => '*.cdna.*.fa.gz' },
        { pattern => '*.pep.*.fa.gz' },
        { pattern => '*.ncrna.*.fa.gz' },
      ],
      -flow_into  => {
        2 =>  { Index => { file => '#file#'} }
      },
    },
    
    {
      -logic_name => 'Index',
      -module     => 'Bio::EnsEMBL::Pipeline::FASTA::NcbiBlastReIndexer',
      -analysis_capacity => 10,
      -rc_name => 'dump',
      -parameters => {
        program => $self->o('ncbiblast_exe')
      }
    },
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
    %{ $self->SUPER::pipeline_wide_parameters() },
    base_path => $self->o('base_path'),
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
    %{$self->SUPER::resource_classes()},
    dump => { 'LSF' => '-q normal -M4000000 -R"select[mem>4000] rusage[mem=4000]"'},
  }
}

1;
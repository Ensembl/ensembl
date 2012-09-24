package Bio::EnsEMBL::Pipeline::PipeConfig::FASTA_conf;

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
        
        #'registry' => 'Reg.pm', # default option to refer to Reg.pm, should be full path
        #'base_path' => '', #where do you want your files
        
        ### Optional overrides        
        ftp_dir => '',

        species => [],
        
        dump_types => [],
        
        db_types => [],
        
        force_species => [],
        
        process_logic_names => [],
        
        skip_logic_names => [],
        
        release => software_version(),
        
        previous_release => (software_version() - 1),
        
        ### SCP code
        
        blast_servers => [],
        blast_genomic_dir => '',
        blast_genes_dir => '',
        
        scp_user => $self->o('ENV', 'USER'),
        scp_identity => '',
        no_scp => 0,
        
        ### Defaults 
        
        pipeline_name => 'fasta_dump_'.$self->o('release'),
        
        wublast_exe => 'xdformat',
        blat_exe => 'faToTwoBit',
        port_offset => 30000,
        
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
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::ReuseSpeciesFactory',
        -parameters => {
          species => $self->o('species'),
          sequence_type_list => $self->o('dump_types'),
          ftp_dir => $self->o('ftp_dir'),
          force_species => $self->o('force_species'),
        },
        -input_ids  => [ {} ],
        -flow_into  => {
          1 => 'Notify',
          2 => 'DumpDNA',
          3 => 'DumpGenes',
          4 => 'CopyDNA',
          5 => 'ChecksumGeneratorFactory'
        },
      },
      
      ######### DUMPING DATA
      
      {
        -logic_name => 'DumpDNA',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::DumpFile',
        -parameters => {
          process_logic_names => $self->o('process_logic_names'),
          skip_logic_names => $self->o('skip_logic_names'),
        },
        -can_be_empty => 1,
        -flow_into  => {
          1 => 'ConcatFiles'
        },
        -can_be_empty     => 1,
        -max_retry_count  => 1,
        -hive_capacity    => 10,
        -rc_name          => 'dump',
      },
      
      {
        -logic_name => 'DumpGenes',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::DumpFile',
        -flow_into  => {
          2 => ['BlastPepIndex'],
          3 => ['BlastGeneIndex']
        },
        -max_retry_count  => 1,
        -hive_capacity    => 10,
        -can_be_empty     => 1,
        -rc_name          => 'dump',
        -wait_for         => 'DumpDNA' #block until DNA is done
      },
      
      {
        -logic_name => 'ConcatFiles',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::ConcatFiles',
        -can_be_empty => 1,
        -max_retry_count => 5,
        -flow_into  => {
          1 => [qw/BlastDNAIndex BlatDNAIndex BlatSmDNAIndex PrimaryAssembly/]
        },
      },
      
      {
        -logic_name       => 'PrimaryAssembly',
        -module           => 'Bio::EnsEMBL::Pipeline::FASTA::CreatePrimaryAssembly',
        -can_be_empty     => 1,
        -max_retry_count  => 5,
      },
      
      ######## COPY DATA
      
      {
        -logic_name => 'CopyDNA',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::CopyDNA',
        -can_be_empty => 1,
        -hive_capacity => 5,
        -parameters => {
          ftp_dir => $self->o('ftp_dir')
        },
      },
      
      ######## INDEXING
      
      {
        -logic_name => 'BlastDNAIndex',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::WuBlastIndexer',
        -parameters => {
          molecule => 'dna', type => 'genomic', program => $self->o('wublast_exe')
        },
        -hive_capacity => 10,
        -can_be_empty => 1,
        -rc_name => 'indexing',
      },
      
      {
        -logic_name => 'BlastPepIndex',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::WuBlastIndexer',
        -parameters => {
          molecule => 'pep', type => 'genes', program => $self->o('wublast_exe')
        },
        -hive_capacity => 5,
        -can_be_empty => 1,
        -flow_into => {
          1 => [qw/SCPBlast/],
        },
      },
      
      {
        -logic_name => 'BlastGeneIndex',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::WuBlastIndexer',
        -parameters => {
          molecule => 'dna', type => 'genes', program => $self->o('wublast_exe')
        },
        -hive_capacity => 5,
        -can_be_empty => 1,
        -flow_into => {
          1 => [qw/SCPBlast/],
        },
      },
      
      {
        -logic_name => 'BlatDNAIndex',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::BlatIndexer',
        -parameters => {
          port_offset => $self->o('port_offset'), 
          program => $self->o('blat_exe'),
          'index' => 'dna' 
        },
        -can_be_empty => 1,
        -hive_capacity => 5,
        -rc_name => 'indexing',
      },
      
      {
        -logic_name => 'BlatSmDNAIndex',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::BlatIndexer',
        -parameters => {
          port_offset => $self->o('port_offset'), 
          program => $self->o('blat_exe'),
          'index' => 'dna_sm' 
        },
        -can_be_empty => 1,
        -hive_capacity => 5,
        -rc_name => 'indexing',
      },
      
      ######## COPYING
      {
        -logic_name => 'SCPBlast',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::SCPBlast',
        -parameters => {
          target_servers => $self->o('blast_servers'),
          genomic_dir => $self->o('blast_genomic_dir'),
          genes_dir => $self->o('blast_genes_dir'),
          
          scp_user => $self->o('scp_user'),
          scp_identity => $self->o('scp_identity'),
          
          no_scp => $self->o('no_scp'),
        },
        -hive_capacity => 3,
        -can_be_empty => 1,
        -wait_for => [qw/DumpDNA DumpGenes PrimaryAssembly BlastDNAIndex BlastGeneIndex BlastPepIndex/]
      },
      
      ####### CHECKSUMMING
      
      {
        -logic_name => 'ChecksumGeneratorFactory',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::FindDirs',
        -parameters => {
          column_names => [qw/dir/],
          input_id => { 'dir' => '#dir#' },
          fan_branch_code => 2,
        },
        -wait_for   => [qw/DumpDNA DumpGenes PrimaryAssembly BlastDNAIndex BlastGeneIndex BlastPepIndex/],
        -flow_into  => { 2 => ['ChecksumGenerator'] } 
      },
      
      {
        -logic_name => 'ChecksumGenerator',
        -module     => 'Bio::EnsEMBL::Pipeline::ChecksumGenerator',
        -hive_capacity => 10,
      },
      
      ####### NOTIFICATION
      
      {
        -logic_name => 'Notify',
        -module     => 'Bio::EnsEMBL::Pipeline::FASTA::EmailSummary',
        -parameters => {
          email   => $self->o('email'),
          subject => $self->o('pipeline_name').' has finished',
        },
        -wait_for   => ['SCPBlast', 'ChecksumGenerator'],
      }
    
    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    
    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
        base_path => $self->o('base_path'), 
        db_types => $self->o('db_types'),
        release => $self->o('release'),
        previous_release => $self->o('previous_release'),
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
      'dump'      => { LSF => '-q long -M1000000 -R"select[mem>1000] rusage[mem=1000]"' },
      'indexing'  => { LSF => '-q normal -M2000000 -R"select[mem>2000] rusage[mem=2000]"' },
    }
}

1;

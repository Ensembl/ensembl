package Bio::EnsEMBL::Pipeline::FASTA::EmailSummary;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Hive::RunnableDB::NotifyByEmail Bio::EnsEMBL::Pipeline::Base/;
use Bio::EnsEMBL::Hive::Utils qw/destringify/;

sub fetch_input {
  my ($self) = @_;
  
  $self->assert_executable('sendmail');
  
  my $dump_dna = $self->jobs('DumpDNA');
  my $copy_dna = $self->jobs('CopyDNA');
  my $dump_genes = $self->jobs('DumpGenes');
  my $blast_dna = $self->jobs('BlastDNAIndex');
  my $blast_gene = $self->jobs('BlastGeneIndex');
  my $blast_pep = $self->jobs('BlastPepIndex');
  my $blat = $self->jobs('BlatDNAIndex', 100);
  my $blat_sm = $self->jobs('BlatSmDNAIndex', 100);
    
  my @args = (
    $dump_dna->{successful_jobs},
    $dump_dna->{failed_jobs},
    $copy_dna->{successful_jobs},
    $copy_dna->{failed_jobs},
    $dump_genes->{successful_jobs},
    $dump_genes->{failed_jobs},
    $blast_dna->{successful_jobs},
    $blast_dna->{failed_jobs},
    $blast_gene->{successful_jobs},
    $blast_gene->{failed_jobs},
    $blast_pep->{successful_jobs},
    $blast_pep->{failed_jobs},
    $blat->{successful_jobs},
    $blat->{failed_jobs},
    $blat_sm->{successful_jobs},
    $blat_sm->{failed_jobs},
    $self->failed(),
    $self->summary($dump_dna),
    $self->summary($copy_dna),
    $self->summary($dump_genes),
    $self->summary($blast_dna),
    $self->summary($blast_gene),
    $self->summary($blast_pep),
    $self->summary($blat),
    $self->summary($blat_sm),
  );
  
  my $msg = sprintf(<<'MSG', @args);
Your FASTA Pipeline has finished. We have:

  * %d species with new DNA Dumped (%d failed)
  * %d species with copied DNA (%d failed)
  * %d species with genes dumped (%d failed)
  * %d species with BLAST DNA indexes generated (%d failed)
  * %d species with BLAST GENE indexes generated (%d failed)
  * %d species with BLAST PEPTIDE indexes generated (%d failed)
  * %d species with BLAT DNA generated (%d failed)
  * %d species with BLAT soft-masked DNA generated (%d failed)

%s

===============================================================================

Full breakdown follows ...

%s

%s

%s

%s

%s

%s

%s

%s

MSG
  $self->param('text', $msg);
  return;
}

sub jobs {
  my ($self, $logic_name, $minimum_runtime) = @_;
  my $aa = $self->db->get_AnalysisAdaptor();
  my $aja = $self->db->get_AnalysisJobAdaptor();
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  my $id = $analysis->dbID();
  my @jobs = @{$aja->generic_fetch("j.analysis_id =$id")};
  $_->{input} = destringify($_->input_id()) for @jobs;
  @jobs = 
    sort { $a->{input}->{species} cmp $b->{input}->{species} }
    grep { 
      if($minimum_runtime) {
        if($minimum_runtime > $_->runtime_msec()) {
          1;
        }
        else {
          0;
        }
      }
      else {
        1;
      }
    }
    @jobs;
  my %passed_species = map { $_->{input}->{species}, 1 } grep { $_->status() eq 'DONE' } @jobs;
  my %failed_species = map { $_->{input}->{species}, 1 } grep { $_->status() eq 'FAILED' } @jobs;
  return {
    analysis => $analysis,
    name => $logic_name,
    jobs => \@jobs,
    successful_jobs => scalar(keys %passed_species),
    failed_jobs => scalar(keys %failed_species),
  };
}


sub failed {
  my ($self) = @_;
  my $failed = $self->db()->get_AnalysisJobAdaptor()->fetch_all_failed_jobs();
  if(! @{$failed}) {
    return 'No jobs failed. Congratulations!';
  }
  my $output = <<'MSG';
The following jobs have failed during this run. Please check your hive's error msg table for the following jobs:

MSG
  foreach my $job (@{$failed}) {
    my $analysis = $self->db()->get_AnalysisAdaptor()->fetch_by_dbID($job->analysis_id());
    my $line = sprintf(q{  * job_id=%d %s(%5d) input_id='%s'}, $job->dbID(), $analysis->logic_name(), $analysis->dbID(), $job->input_id());
    $output .= $line;
    $output .= "\n";
  }
  return $output;
}

my $sorter = sub {
  my $status_to_int = sub {
    my ($v) = @_;
    return ($v->status() eq 'FAILED') ? 0 : 1;
  };
  my $status_sort = $status_to_int->($a) <=> $status_to_int->($b);
  return $status_sort if $status_sort != 0;
  return $a->{input}->{species} cmp $b->{input}->{species};
};

sub summary {
  my ($self, $data) = @_;
  my $name = $data->{name};
  my $underline = '~'x(length($name));
  my $output = "$name\n$underline\n\n";
  my @jobs = @{$data->{jobs}};
  if(@jobs) {
    foreach my $job (sort $sorter @{$data->{jobs}}) {
      my $species = $job->{input}->{species};
      $output .= sprintf("  * %s - job_id=%d %s\n", $species, $job->dbID(), $job->status());
    }
  }
  else {
    $output .= "No jobs run for this analysis\n";
  }
  $output .= "\n";
  return $output;
}

1;
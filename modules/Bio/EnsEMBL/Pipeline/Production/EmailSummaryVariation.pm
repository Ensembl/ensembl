package Bio::EnsEMBL::Pipeline::Production::EmailSummaryVariation;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Hive::RunnableDB::NotifyByEmail Bio::EnsEMBL::Pipeline::Base/;
use Bio::EnsEMBL::Hive::Utils qw/destringify/;

sub fetch_input {
  my ($self) = @_;
  
  $self->assert_executable('sendmail');
  
  my $snp_density = $self->jobs('SnpDensity');
  my $snp_count = $self->jobs('SnpCount');
  my $nonsense = $self->jobs('NonSense');
    
  my @args = (
    $snp_density->{successful_jobs},
    $snp_density->{failed_jobs},
    $snp_count->{successful_jobs},
    $snp_count->{failed_jobs},
    $nonsense->{successful_jobs},
    $nonsense->{failed_jobs},
    $self->failed(),
    $self->summary($snp_density),
    $self->summary($snp_count),
    $self->summary($nonsense),
  );
  
  my $msg = sprintf(<<'MSG', @args);
Your FASTA Pipeline has finished. We have:

  * %d species with snp density (%d failed)
  * %d species with snp count (%d failed)
  * %d species with nonsense attributes (%d failed)

%s

===============================================================================

Full breakdown follows ...

%s

%s

%s

MSG
  $self->param('text', $msg);
  return;
}

sub jobs {
  my ($self, $logic_name) = @_;
  my $aa = $self->db->get_AnalysisAdaptor();
  my $aja = $self->db->get_AnalysisJobAdaptor();
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  my @jobs;
  if (!$analysis) {
    return {
      name => $logic_name,
      successful_jobs => 0,
      failed_jobs => 0,
      jobs => \@jobs,
    };
  }
  my $id = $analysis->dbID();
  @jobs = @{$aja->generic_fetch("j.analysis_id =$id")};
  $_->{input} = destringify($_->input_id()) for @jobs;
  @jobs = sort { $a->{input}->{species} cmp $b->{input}->{species} } @jobs;
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

# Base class that all other mapping methods inherit from

package XrefMapper::Methods::ExonerateBasic;

use strict;

use File::Basename;
use IPC::Open3;

# Path to exonerate executable
my $exonerate_path = "/usr/local/ensembl/bin/exonerate-0.8.3";


sub new {

  my($class) = @_;

  my $self ={};
  bless $self,$class;

  return $self;

}

=head2 run

  Arg[1]     : Query filename to pass to exonerate; this should be the XREF file.
  Arg[2]     : Target filename to pass to exonerate; this should be the ENSEMBL file.
  Example    : none
  Description: Run exonerate with default options.
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub run() {

  my ($self, $query, $target) = @_;

  my $name = $self->submit_exonerate($query, $target, options());

  return $name;

}

=head2 options

  Args       : none
  Example    : none
  Description: Options to pass to exonerate. Override as appropriate.
  Returntype : List of strings.
  Exceptions : none
  Caller     : general

=cut

sub options() {

  return ();

}

=head2 submit_exonerate

  Args       : none
  Example    : none
  Description: Submit a *single* mapping job array to exonerate.
  Returntype : The name of the job array submitted to LSF. Can be used in depend job.
  Exceptions : none
  Caller     : general

=cut

sub submit_exonerate {

  my ($self, $query, $target, @options) = @_;

  my $num_jobs = calculate_num_jobs($query);

  # TODO - get root_dir from config
  my $root_dir = "/nfs/acari/gp1/work/ensembl/misc-scripts/xref_mapping";

  my $options_str = join(" ", @options);

  my $unique_name = $self->get_class_name() . "_" . time();

  my $prefix = $root_dir . "/" . basename($query);
  $prefix =~ s/\.\w+$//;

  my $output = $self->get_class_name() . "_\$LSB_JOBINDEX.map";

  my @main_bsub = ( 'bsub', '-J', $unique_name . "[1-$num_jobs]", '-o', "$prefix.%I.out", '-e', "$prefix.%I.err");

  #print "bsub command: " . join(" ", @main_bsub) . "\n\n";

  # Use IPC::Open3 to open the process so that we can read and write from/to its stdout/stderr
  my ($wtr, $rtr, $etr, $pid);
  $pid = open3($wtr, $rtr, $etr, @main_bsub);

  # Create actual execute script to be executed with LSF, and write to pipe
  my $main_job = <<EOF;
. /usr/local/lsf/conf/profile.lsf

cd /tmp

rm -f /tmp/\$LSB_JOBINDEX.$query /tmp/\$LSB_JOBINDEX.$target /tmp/$output

lsrcp ecs1a:$root_dir/$target /tmp/\$LSB_JOBINDEX.$target
lsrcp ecs1a:$root_dir/$query  /tmp/\$LSB_JOBINDEX.$query

$exonerate_path /tmp/\$LSB_JOBINDEX.$query /tmp/\$LSB_JOBINDEX.$target --querychunkid \$LSB_JOBINDEX --querychunktotal $num_jobs --showvulgar false --showalignment FALSE --ryo "xref:%qi:%ti:%qab:%qae:%tab:%tae:%C:%s\n" $options_str | grep '^xref' > /tmp/$output

lsrcp /tmp/$output ecs1a:$root_dir/$output

rm -f /tmp/\$LSB_JOBINDEX.$query /tmp/\$LSB_JOBINDEX.$target /tmp/$output
EOF

  # TODO make sure query/target are the right way round

  print $wtr $main_job;

  #print "Job:\n" . $main_job . "\n\n";

  close($wtr);

  # Wait until bsub has actually run - will print to its stdout ($rtr) and then close it
  my $jobid;
  while (<$rtr>) {
    if (/Job <([0-9]+)> is/) {
      $jobid = $1;
      print "LSF job ID for main mapping job: $jobid \n"
    }
  }
  if (!defined($jobid)) {
    print STDERR "Error: could not get LSF job ID for mapping job\n";
  }

  return $unique_name;

}

=head2 calculate_num_jobs

  Args       : Query file name
  Example    : none
  Description: Calculate the number of LSF jobs to submit based on the size of the query file.
  Returntype : The number of jobs.
  Exceptions : none
  Caller     : general

=cut

sub calculate_num_jobs {

  my $query = shift;

  my $bytes_per_job = 250000;

  my $size = (stat $query)[7];

  return int($size/$bytes_per_job);

}

# Get class name from fully-qualified object name
# e.g. return ExonerateBasic from XrefMapper::Methods::ExonerateBasic=HASH(Ox12113c0)

sub get_class_name() {

  my $self = shift;

  $self =~ s/=.*$//;

  my @bits = split(/::/, $self);

  return @bits[$#bits];

}

1;


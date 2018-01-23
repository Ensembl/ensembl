=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

package XrefMapper::Methods::PBSExonerateGappedBest1;

use strict;
#use File::Basename;
#use IPC::Open3;
use Sys::Hostname;

use XrefMapper::Methods::ExonerateBasic;

use vars '@ISA';

@ISA = qw{XrefMapper::Methods::ExonerateBasic};

my $exonerate_path = "exonerate";

sub options {

  return ('--model', 'affine:local', '--subopt', 'no', '--bestn', '1');

}

sub query_identity_threshold {

  return 90;

}

sub target_identity_threshold {

  return 90;

}

sub submit_exonerate{
  my ($self, $query, $target, $root_dir, @options) = @_;

  # Exonerate can run individual chunks of a large job.
  # Determine the number of chunks that will be used, and add to total
  my $num_jobs 
      = XrefMapper::Methods::ExonerateBasic::calculate_num_jobs($query);
  $self->jobcount($self->jobcount()+$num_jobs);

  $num_jobs || return; # Check we have jobs to run

  # Build a template representing the exonerate command
  my $options_str = join( ' ', @options );
  my $shell_command = <<"EOF";
#!/bin/sh
$exonerate_path \\
--target $target \\
--query $query \\
--querychunktotal $num_jobs \\
--querychunkid %d \\
--showvulgar FALSE \\
--showalignment FALSE \\
--ryo "xref:%%qi:%%ti:%%ei:%%ql:%%tl:%%qab:%%qae:%%tab:%%tae:%%C:%%s\\n" \\
$options_str \\
| grep '^xref'
exit
EOF

  #warn( sprintf $shell_command, 1 );

  # Determine the output filename
  my $host = Sys::Hostname::hostname;
  my $query_alphabet = $query =~ /peptide.fasta$/ ? 'peptide' : 'dna';
  my $outfile_root = sprintf( "%s:%s/%s_%s",
                              $host,
                              $root_dir,
                              $self->get_class_name(), 
                              $query_alphabet );
  
  # Run PBS qsub for each job chunk.
  # First set a job that depends on the completion of all other jobs
  my $jobid = $self->submit_qsub( q( echo 'sleep 1' ), 
                                  [ -W => "depend=on:$num_jobs" ] );

  # Run each job chunk
  my $jobname = join '.', 'X', ( $host =~ /^(\w+)/ ), $$; 
  foreach( my $i=1; $i<=$num_jobs; $i++ ){

    my $outfile = sprintf( '%s_%s.map', $outfile_root, $i );
    my $errfile = sprintf( '%s_%s.err', $outfile_root, $i );

    my $chunkid = $self->submit_qsub
        ( sprintf( $shell_command, $i ),
          [ '-N' => $jobname, 
            '-W' => "depend=beforeany:$jobid", 
            '-o' => $outfile, 
            '-e' => $errfile ] );
  }

  # BasicMapper calls a method called 'submit_depend_job' that waits
  # till all the lsf jobs have completed. We're not using LSF, so this
  # approach does not work. Solution: wait here until all jobs are done.
  # Not as efficient for a large cluster, but we may have enough work to
  # do to saturate a small one. TODO: abstract 'submit_depend_job' to the 
  # Method class instances.

  $self->global_jobid( $jobid );
  #$self->submit_depend_job( $jobid );

  return $jobid;
}

#----------------------------------------------------------------------
# Get/set the 'global' job ID that groups all chunked jobs in the method
sub global_jobid{
  my $key = '_global_jobid';
  my $self = shift;
  if( @_ ){ $self->{$key} = shift }
  return $self->{$key};
}

#----------------------------------------------------------------------
# Waits in-process untill the specified jobid has completed
sub submit_depend_job{
  my $self = shift;

  $self->jobcount || return; # Check that we are running jobs
  my $jobid = $self->global_jobid || die( "global_jobid unset!" );

  my $depend = join( ':', 'afterany', $jobid );
  
  my $id = $self->submit_qsub( q( echo 'sleep 1' ), 
                               [ '-W' => 'block=true',
                                 '-W' => "depend=$depend" ] );
  return $jobid;
}

#----------------------------------------------------------------------
# A wrapper for submitting PBS qsub jobs;
# First the READER (this process) opens a pipe (-|) on the WRITER.
# The WRITER then opens a pipe (|-) on the RUNNER.
# The RUNNER then execs qsub with command line options,
# the WRITER writes the script to the RUNNER, 
# and any output is collected by the READER.
# Returns the qsub job ID.
sub submit_qsub{
  my $self = shift;
  my $script = shift || die( "Need a script to submit to qsub!" );
  my @qsub_args = @{ shift || [] };

  local *QSUB;
  local *QSUB_READER;
  
  my $jobid;
  
  my $writer_pid;
  if( $writer_pid = open( QSUB_READER, '-|' ) ){
    # READER - Reads stdout from RUNNER process
    while( <QSUB_READER> ){
      if( /^(\d+)/ ){
        $jobid = $1;
        print( "Job ID $1 submitted to qsub\n" );
      }
    }
    close( QSUB_READER );
  }
  else{
    unless( defined($writer_pid) ){ die("Could not fork : $!" ) }
    my $runner_pid;
    if( $runner_pid = open(QSUB, '|-') ){
      # WRITER - Writes command to RUNNER process
      print QSUB $script;
      close QSUB;
      unless( $? == 0 ){ die( "qsub failed; non-zero status" ) }
      exit(0);
    }
    else{
      # RUNNER - Runs the command; STDIN,STDOUT attached to WRITER,READER.
      unless( defined($runner_pid) ){ die("Could not fork : $!" ) }
      #warn join( " ", 'qsub', @qsub_args );
      exec( 'qsub', @qsub_args );
      die("Could not exec qsub : $!");
    }
  }
  return $jobid;
}


1;

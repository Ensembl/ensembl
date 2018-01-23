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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Utils::ConversionSupport - Utility module for Vega release and
schema conversion scripts

=head1 SYNOPSIS

  my $serverroot = '/path/to/ensembl';
  my $suport = new Bio::EnsEMBL::Utils::ConversionSupport($serverroot);

  # parse common options
  $support->parse_common_options;

  # parse extra options for your script
  $support->parse_extra_options( 'string_opt=s', 'numeric_opt=n' );

  # ask user if he wants to run script with these parameters
  $support->confirm_params;

  # see individual method documentation for more stuff

=head1 DESCRIPTION

This module is a collection of common methods and provides helper
functions for the Vega release and schema conversion scripts. Amongst
others, it reads options from a config file, parses commandline options
and does logging.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Logger;

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin $Script);
use POSIX qw(strftime);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::ScriptUtils qw(parse_bytes);
use Carp;

my %level_defs = (
  'error'     => 1,
  'warn'      => 2,
  'warning'   => 2,
  'info'      => 3,
  'debug'     => 4,
  'verbose'   => 4,
);

my @reverse_level_defs = (undef, qw(error warning info debug));

=head2 new

  Arg[1]      : String $serverroot - root directory of your ensembl sandbox
  Example     : my $support = new Bio::EnsEMBL::Utils::ConversionSupport(
                                        '/path/to/ensembl');
  Description : constructor
  Return type : Bio::EnsEMBL::Utils::ConversionSupport object
  Exceptions  : thrown on invalid loglevel
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($logfile, $logauto, $logautobase, $logautoid, $logpath, $logappend,
    $loglevel, $is_component) = rearrange(
      ['LOGFILE', 'LOGAUTO', 'LOGAUTOBASE', 'LOGAUTOID', 'LOGPATH', 'LOGAPPEND',
      'LOGLEVEL', 'IS_COMPONENT'], @_);
  
  my $self = { '_warnings'     => 0, };
  bless ($self, $class);

  # initialise
  $self->logfile($logfile);
  $self->logpath($logpath);
  $self->logappend($logappend);
  $self->is_component($is_component);

  # automatic logfile creation
  $self->logauto($logauto);
  $logautoid ||= strftime("%Y%m%d-%H%M%S", localtime);
  $self->log_auto_id($logautoid);
  $self->create_auto_logfile($logautobase);

  $loglevel ||= 'info';
  if ($loglevel =~ /^\d+$/ and $loglevel > 0 and $loglevel < 5) {
    $self->{'loglevel'} = $loglevel;
  } elsif ($level_defs{lc($loglevel)}) {
    $self->{'loglevel'} = $level_defs{lc($loglevel)};
  } else {
    throw('Unknown loglevel: $loglevel.');
  }
  
  return $self;
}


=head2 log_generic

  Arg[1]      : String $txt - the text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle;
                $support->log('Log foo.\n', 1);
  Description : Logs a message to the filehandle initialised by calling
                $self->log_filehandle(). You can supply an indentation level
                to get nice hierarchical log messages.
  Return type : true on success
  Exceptions  : thrown when no filehandle can be obtained
  Caller      : general

=cut

sub log_generic {
  my ($self, $txt, $indent, $stamped) = @_;

  $indent ||= 0;
  my $fh = $self->log_filehandle;

  # append timestamp and memory usage to log text if requested
  if ($stamped) {
    $txt =~ s/^(\n*)(.*)(\n*)$/$2/;
    $txt = sprintf("%-60s%20s", $txt, $self->time_and_mem);
    $txt = $1.$txt.$3;
  }
  
  # strip off leading linebreaks so that indenting doesn't break
  $txt =~ s/^(\n*)//;
  
  # indent
  $txt = $1."  "x$indent . $txt;
  
  print $fh "$txt";
  
  return(1);
}


=head2 error

  Arg[1]      : String $txt - the error text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle;
                $support->log_error('Log foo.\n', 1);
  Description : Logs a message via $self->log and exits the script.
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub error {
  my ($self, $txt, $indent, $stamped) = @_;

  return(0) unless ($self->{'loglevel'} >= 1);
  
  $txt = "ERROR: ".$txt;
  $self->log_generic($txt, $indent, $stamped);
  
  $self->log_generic("\nExiting prematurely.\n\n");
  $self->log_generic("Runtime: ".$self->runtime." ".$self->date_and_mem."\n\n");
  
  confess($txt);
}


=head2 warning

  Arg[1]      : String $txt - the warning text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle;
                $support->log_warning('Log foo.\n', 1);
  Description : Logs a message via $self->log and increases the warning counter.
  Return type : true on success
  Exceptions  : none
  Caller      : general

=cut

sub warning {
  my ($self, $txt, $indent, $stamped) = @_;
  
  return(0) unless ($self->{'loglevel'} >= 2);
  
  $txt = "WARNING: " . $txt;
  $self->log_generic($txt, $indent, $stamped);
  
  $self->{'_warnings'}++;
  
  return(1);
}


sub info {
  my ($self, $txt, $indent, $stamped) = @_;

  return(0) unless ($self->{'loglevel'} >= 3);

  $self->log_generic($txt, $indent, $stamped);
  return(1);
}


=head2 debug

  Arg[1]      : String $txt - the warning text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle;
                $support->log_verbose('Log this verbose message.\n', 1);
  Description : Logs a message via $self->log if --verbose option was used
  Return type : TRUE on success, FALSE if not verbose
  Exceptions  : none
  Caller      : general

=cut

sub debug {
  my ($self, $txt, $indent, $stamped) = @_;

  return(0) unless ($self->{'loglevel'} >= 4);

  $self->log_generic($txt, $indent, $stamped);
  return(1);
}


sub log_progress {
  my $self = shift;
  my $name = shift;
  my $curr = shift;
  my $indent = shift;
  
  throw("You must provide a name and the current value for your progress bar")
    unless ($name and $curr);

  # return if we haven't reached the next increment
  return if ($curr < int($self->{'_progress'}->{$name}->{'next'}));

  my $index = $self->{'_progress'}->{$name}->{'index'};
  my $num_bins = $self->{'_progress'}->{$name}->{'numbins'};
  my $percent = $index/$num_bins*100;

  my $log_str;
  $log_str .= '      'x$indent if ($index == 0);
  $log_str .= "\b"x4;
  $log_str .= sprintf("%3s%%", $percent);
  $log_str .= "\n" if ($curr == $self->{'_progress'}->{$name}->{'max_val'});

  $self->info($log_str);

  # increment counters
  $self->{'_progress'}->{$name}->{'index'}++;
  $self->{'_progress'}->{$name}->{'next'} += $self->{'_progress'}->{$name}->{'binsize'};
}


sub log_progressbar {
  my $self = shift;
  my $name = shift;
  my $curr = shift;
  my $indent = shift;
  
  throw("You must provide a name and the current value for your progress bar")
    unless ($name and $curr);

  # return if we haven't reached the next increment
  return if ($curr < int($self->{'_progress'}->{$name}->{'next'}));

  my $index = $self->{'_progress'}->{$name}->{'index'};
  my $num_bins = $self->{'_progress'}->{$name}->{'numbins'};
  my $percent = $index/$num_bins*100;

  my $log_str = "\r".('  'x$indent)."[".('='x$index).(' 'x($num_bins-$index))."] ${percent}\%";
  $log_str .= "\n" if ($curr == $self->{'_progress'}->{$name}->{'max_val'});

  $self->info($log_str);

  # increment counters
  $self->{'_progress'}->{$name}->{'index'}++;
  $self->{'_progress'}->{$name}->{'next'} += $self->{'_progress'}->{$name}->{'binsize'};
}


sub init_progress {
  my $self = shift;
  my $max = shift;
  my $num_bins = shift || 50;

  throw("You must provide the maximum value for your progress bar")
    unless (defined($max));

  # auto-generate a unique name for your progressbar
  my $name = time . '_' . int(rand(1000));

  # calculate bin size; we will use 50 bins (2% increments)
  my $binsize = $max/$num_bins;

  $self->{'_progress'}->{$name}->{'max_val'} = $max;
  $self->{'_progress'}->{$name}->{'binsize'} = $binsize;
  $self->{'_progress'}->{$name}->{'numbins'} = $num_bins;
  $self->{'_progress'}->{$name}->{'next'} = 0;
  $self->{'_progress'}->{$name}->{'index'} = 0;

  return $name;
}


=head2 log_filehandle

  Arg[1]      : (optional) String $mode - file access mode
  Example     : my $log = $support->log_filehandle;
                # print to the filehandle
                print $log 'Lets start logging...\n';
                # log via the wrapper $self->log()
                $support->log('Another log message.\n');
  Description : Returns a filehandle for logging (STDERR by default, logfile if
                set from config or commandline). You can use the filehandle
                directly to print to, or use the smart wrapper $self->log().
                Logging mode (truncate or append) can be set by passing the
                mode as an argument to log_filehandle(), or with the
                --logappend commandline option (default: truncate)
  Return type : Filehandle - the filehandle to log to
  Exceptions  : thrown if logfile can't be opened
  Caller      : general

=cut

sub log_filehandle {
  my ($self, $mode) = @_;
  
  unless ($self->{'_log_filehandle'}) {
    $mode ||= '>';
    $mode = '>>' if ($self->logappend);
    
    my $fh = \*STDERR;
    
    if (my $logfile = $self->logfile) {
      if (my $logpath = $self->logpath) {
        unless (-e $logpath) {
          system("mkdir -p $logpath") == 0 or
            throw("Can't create log dir $logpath: $!\n");
        }
        
        $logfile = "$logpath/".$self->logfile;
      }
      
      open($fh, "$mode", $logfile) or
        throw("Unable to open $logfile for writing: $!");
    }

    $self->{'_log_filehandle'} = $fh;
  }

  return $self->{'_log_filehandle'};
}


=head2 extract_log_identifier 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub extract_log_identifier {
  my $self = shift;

  if (my $logfile = $self->logfile) {
    $logfile =~ /.+\.([^\.]+)\.log/;
    return $1;
  } else {
    return undef;
  }
}


=head2 init_log

  Example     : $support->init_log;
  Description : Opens a filehandle to the logfile and prints some header
                information to this file. This includes script name, date, user
                running the script and parameters the script will be running
                with.
  Return type : Filehandle - the log filehandle
  Exceptions  : none
  Caller      : general

=cut

sub init_log {
  my $self = shift;
  my $params = shift;

  # get a log filehandle
  my $log = $self->log_filehandle;

  # remember start time
  $self->{'_start_time'} = time;

  # don't log parameters if this script is run by another one
  if ($self->logauto or ! $self->is_component) {
    # print script name, date, user who is running it
    my $hostname = `hostname`;
    chomp $hostname;
    my $script = "$hostname:$Bin/$Script";
    my $user = `whoami`;
    chomp $user;
    $self->info("Script: $script\nDate: ".$self->date."\nUser: $user\n");

    # print parameters the script is running with
    if ($params) {
      $self->info("Parameters:\n\n");
      $self->info($params);
    }
  }

  return $log;
}


=head2 finish_log

  Example     : $support->finish_log;
  Description : Writes footer information to a logfile. This includes the
                number of logged warnings, timestamp and memory footprint.
  Return type : TRUE on success
  Exceptions  : none
  Caller      : general

=cut

sub finish_log {
  my $self = shift;
  
  $self->info("\nAll done for $Script.\n");
  $self->info($self->warning_count." warnings. ");
  $self->info("Runtime: ".$self->runtime." ".$self->date_and_mem."\n\n");
  
  return(1);
}


sub runtime {
  my $self = shift;

  my $runtime = "n/a";

  if ($self->{'_start_time'}) {
    my $diff = time - $self->{'_start_time'};
    my $sec = $diff % 60;
    $diff = ($diff - $sec) / 60;
    my $min = $diff % 60;
    my $hours = ($diff - $min) / 60;
    
    $runtime = "${hours}h ${min}min ${sec}sec";
  }

  return $runtime;
}


=head2 date_and_mem

  Example     : print LOG "Time, memory usage: ".$support->date_and_mem."\n";
  Description : Prints a timestamp and the memory usage of your script.
  Return type : String - timestamp and memory usage
  Exceptions  : none
  Caller      : general

=cut

sub date_and_mem {
  my $date = strftime "%Y-%m-%d %T", localtime;
  my $mem = `ps -p $$ -o vsz |tail -1`;
  chomp $mem;
  $mem = parse_bytes($mem*1000);
  return "[$date, mem $mem]";
}


sub time_and_mem {
  my $date = strftime "%T", localtime;
  my $mem = `ps -p $$ -o vsz |tail -1`;
  chomp $mem;
  $mem = parse_bytes($mem*1000);
  $mem =~ s/ //;
  return "[$date|$mem]";
}


=head2 date

  Example     : print "Date: " . $support->date . "\n";
  Description : Prints a nicely formatted timestamp (YYYY-DD-MM hh:mm:ss)
  Return type : String - the timestamp
  Exceptions  : none
  Caller      : general

=cut

sub date {
  return strftime "%Y-%m-%d %T", localtime;
}


=head2 mem

  Example     : print "Memory usage: " . $support->mem . "\n";
  Description : Prints the memory used by your script. Not sure about platform
                dependence of this call ...
  Return type : String - memory usage
  Exceptions  : none
  Caller      : general

=cut

sub mem {
  my $mem = `ps -p $$ -o vsz |tail -1`;
  chomp $mem;
  return $mem;
}


=head2 warning_count

  Example     : print LOG "There were ".$support->warnings." warnings.\n";
  Description : Returns the number of warnings encountered while running the
                script (the warning counter is increased by $self->log_warning).
  Return type : Int - number of warnings
  Exceptions  : none
  Caller      : general

=cut

sub warning_count {
  my $self = shift;
  return $self->{'_warnings'};
}


=head2 logfile

  Arg[1]      : 
  Example     : 
  Description : Getter and setter for the logfile
  Return type :
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub logfile {
  my $self = shift;
  $self->{'_logfile'} = shift if (@_);
  return $self->{'_logfile'};
}


=head2 log_auto_id

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub log_auto_id {
  my $self = shift;
  $self->{'_log_auto_id'} = shift if (@_);
  return $self->{'_log_auto_id'};
}


sub logauto {
  my $self = shift;
  $self->{'_log_auto'} = shift if (@_);
  return $self->{'_log_auto'};
}


=head2 create_auto_logfile

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      : At Risk
              : under development

=cut

sub create_auto_logfile {
  my $self = shift;
  my $logautobase = shift;

  # do nothing if automatic logfile generation isn't set
  return unless ($self->logauto);

  # an explicit logfile name overrides LOGAUTO
  return if ($self->logfile);

  # argument check
  unless ($logautobase) {
    throw('Need a base logfile name for auto-generating logfile.');
  }

  # create a logfile name
  $self->logfile("${logautobase}_".$self->log_auto_id.".log");
}


=head2 logpath

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub logpath {
  my $self = shift;
  $self->{'_logpath'} = shift if (@_);
  return $self->{'_logpath'};
}


=head2 logappend

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub logappend {
  my $self = shift;
  $self->{'_logappend'} = shift if (@_);
  return $self->{'_logappend'};
}


=head2 is_component

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub is_component {
  my $self = shift;
  $self->{'_is_component'} = shift if (@_);
  return $self->{'_is_component'};
}


sub loglevel {
  my $self = shift;
  return $reverse_level_defs[$self->{'loglevel'}];
}


#
# deprecated methods (left here for backwards compatibility
#
sub log_error {
  return $_[0]->error(@_);
}

sub log_warning {
  return $_[0]->warning(@_);
}

sub log {
  return $_[0]->info(@_);
}

sub log_verbose {
  return $_[0]->debug(@_);
}

sub log_stamped {
  return $_[0]->log(@_, 1);
}



1;


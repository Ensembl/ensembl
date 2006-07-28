package Bio::EnsEMBL::Utils::Logger;

=head1 NAME

Bio::EnsEMBL::Utils::ConversionSupport - Utility module for Vega release and
schema conversion scripts

=head1 SYNOPSIS

    my $serverroot = '/path/to/ensembl';
    my $suport = new Bio::EnsEMBL::Utils::ConversionSupport($serverroot);

    # parse common options
    $support->parse_common_options;

    # parse extra options for your script
    $support->parse_extra_options('string_opt=s', 'numeric_opt=n');

    # ask user if he wants to run script with these parameters
    $support->confirm_params;

    # see individual method documentation for more stuff

=head1 DESCRIPTION

This module is a collection of common methods and provides helper functions 
for the Vega release and schema conversion scripts. Amongst others, it reads
options from a config file, parses commandline options and does logging.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin $Script);
use POSIX qw(strftime);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 new

  Arg[1]      : String $serverroot - root directory of your ensembl sandbox
  Example     : my $support = new Bio::EnsEMBL::Utils::ConversionSupport(
                                        '/path/to/ensembl');
  Description : constructor
  Return type : Bio::EnsEMBL::Utils::ConversionSupport object
  Exceptions  : thrown if no serverroot is provided
  Caller      : general

=cut

sub new {
  my $class = shift;

  my ($logfile, $logpath, $logappend, $verbose, $is_component) = rearrange(
    ['LOGFILE', 'LOGPATH', 'LOGAPPEND', 'VERBOSE', 'IS_COMPONENT'], @_);
  
  my $self = { '_warnings'     => 0, };
  bless ($self, $class);

  # initialise
  $self->logfile($logfile);
  $self->logpath($logpath);
  $self->logappend($logappend);
  $self->verbose($verbose);
  $self->is_component($is_component);
  
  return $self;
}


=head2 log

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

sub log {
  my ($self, $txt, $indent) = @_;

  $indent ||= 0;
  my $fh = $self->log_filehandle;
  
  # strip off leading linebreaks so that indenting doesn't break
  $txt =~ s/^(\n*)//;
  
  # indent
  $txt = $1."  "x$indent . $txt;
  
  print $fh "$txt";
  
  return(1);
}


=head2 log_warning

  Arg[1]      : String $txt - the warning text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle;
                $support->log_warning('Log foo.\n', 1);
  Description : Logs a message via $self->log and increases the warning counter.
  Return type : true on success
  Exceptions  : none
  Caller      : general

=cut

sub log_warning {
  my ($self, $txt, $indent) = @_;
  
  $txt = "WARNING: " . $txt;
  $self->log($txt, $indent);
  
  $self->{'_warnings'}++;
  
  return(1);
}


=head2 log_error

  Arg[1]      : String $txt - the error text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle;
                $support->log_error('Log foo.\n', 1);
  Description : Logs a message via $self->log and exits the script.
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub log_error {
  my ($self, $txt, $indent) = @_;
  
  $txt = "ERROR: ".$txt;
  $self->log($txt, $indent);
  
  $self->log("\nExiting prematurely.\n\n");
  $self->log("Runtime: ".$self->runtime." ".$self->date_and_mem."\n\n");
  
  exit(1);
}


=head2 log_verbose

  Arg[1]      : String $txt - the warning text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle;
                $support->log_verbose('Log this verbose message.\n', 1);
  Description : Logs a message via $self->log if --verbose option was used
  Return type : TRUE on success, FALSE if not verbose
  Exceptions  : none
  Caller      : general

=cut

sub log_verbose {
  my ($self, $txt, $indent) = @_;

  return(0) unless $self->verbose;

  $self->log($txt, $indent);
  return(1);
}


=head2 log_stamped

  Arg[1]      : String $txt - the warning text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle;
                $support->log_stamped('Log this stamped message.\n', 1);
  Description : Appends timestamp and memory usage to a message and logs it via
                $self->log
  Return type : TRUE on success
  Exceptions  : none
  Caller      : general

=cut

sub log_stamped {
  my ($self, $txt, $indent) = @_;
  
  # append timestamp and memory usage to log text
  $txt =~ s/(\n*)$//;
  $txt .= " ".$self->date_and_mem.$1;
  
  $self->log($txt, $indent);
  return(1);
}


sub log_progress {
  my $self = shift;
  my $max = shift;
  my $curr = shift;
  my $incr = shift || 20;
  my $indent = shift;
  my $show_mem = shift;

  throw("You must provide a maximum and current value to log progress.")
    unless ($max and $curr);

  if (($curr % $incr) == 0 or $curr < 20 or $curr == $max) {
    my $mem;
    $mem = ", mem ".$self->mem if ($show_mem);
    my $log_str = "\r".('  'x$indent)."$curr/$max (".int($curr/$max*100)."\%$mem)";
    $log_str .= "\n" if ($curr == $max);
    
    $self->log($log_str);
  }

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
  my $percent = $index*5;

  my $log_str = "\r".('  'x$indent)."[".('='x$index).(' 'x(20-$index))."] ${percent}\%";
  $log_str .= "\n" if ($curr == $self->{'_progress'}->{$name}->{'max_val'});

  $self->log($log_str);

  # increment counters
  $self->{'_progress'}->{$name}->{'index'}++;
  $self->{'_progress'}->{$name}->{'next'} +=
    $self->{'_progress'}->{$name}->{'binsize'};
}


sub init_progressbar {
  my $self = shift;
  my $name = shift;
  my $max = shift;

  throw("You must provide a name and the maximum value for your progress bar")
    unless ($name and $max);

  # calculate bin size; we will use 20 bins (5% increments)
  my $binsize = $max/20;

  $self->{'_progress'}->{$name}->{'max_val'} = $max;
  $self->{'_progress'}->{$name}->{'binsize'} = $binsize;
  $self->{'_progress'}->{$name}->{'next'} = 0;
  $self->{'_progress'}->{$name}->{'index'} = 0;
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
          system("mkdir $logpath") == 0 or
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


=head2 

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
  unless ($self->is_component) {
    # print script name, date, user who is running it
    my $hostname = `hostname`;
    chomp $hostname;
    my $script = "$hostname:$Bin/$Script";
    my $user = `whoami`;
    chomp $user;
    $self->log("Script: $script\nDate: ".$self->date."\nUser: $user\n");

    # print parameters the script is running with
    if ($params) {
      $self->log("Parameters:\n\n");
      $self->log($params);
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
  
  $self->log("\nAll done for $Script.\n");
  $self->log($self->warnings." warnings. ");
  $self->log("Runtime: ".$self->runtime." ".$self->date_and_mem."\n\n");
  
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
  return "[$date, mem $mem]";
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


=head2 warnings

  Example     : print LOG "There were ".$support->warnings." warnings.\n";
  Description : Returns the number of warnings encountered while running the
                script (the warning counter is increased by $self->log_warning).
  Return type : Int - number of warnings
  Exceptions  : none
  Caller      : general

=cut

sub warnings {
  my $self = shift;
  return $self->{'_warnings'};
}


=head2 

  Arg[1]      : 
  Example     : 
  Description : 
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


=head2 

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


=head2 

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


=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub verbose {
  my $self = shift;
  $self->{'_verbose'} = shift if (@_);
  return $self->{'_verbose'};
}


=head2 

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


1;


package Bio::EnsEMBL::Utils::ConfParser;

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

use Getopt::Long;
use Text::Wrap;
use Cwd qw(abs_path);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(user_proceed);
use Bio::EnsEMBL::Utils::Logger;


=head2 new

  Arg[SERVERROOT] : String $serverroot
                root directory of your ensembl code
  Example     : my $support = new Bio::EnsEMBL::Utils::ConversionSupport(
                                        '/path/to/ensembl');
  Description : constructor
  Return type : Bio::EnsEMBL::Utils::ConversionSupport object
  Exceptions  : thrown if no serverroot is provided
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($serverroot) = rearrange(['SERVERROOT'], @_);

  throw("You must supply a serverroot.") unless ($serverroot);

  my $self = {
      '_serverroot'   => $serverroot,
      '_param'        => { interactive => 1 },
  };

  bless ($self, $class);

  return $self;
}


=head2 parse_common_options

  Example     : $support->parse_common_options;
  Description : This method reads options from a configuration file and parses
                some commandline options that are common to all scripts (like
                db connection settings, help, dry-run). Commandline options
                will override config file settings. 

                All options will be accessible via $self->param('name').
  Return type : true on success 
  Exceptions  : thrown if configuration file can't be opened
  Caller      : general

=cut

sub parse_common_options {
  my $self = shift;

  # read commandline options
  my %h;
  Getopt::Long::Configure("pass_through");
  &GetOptions( \%h,
      'conffile|conf=s',
      'logfile|log=s',
      'logpath=s',
      'logappend|log_append=s',
      'is_component=s',
      'verbose|v=s',
      'interactive|i=s',
      'dry_run|dry|n=s',
      'help|h|?',
      );

  # reads config file
  my $conffile = $h{'conffile'} || $self->param('default_conf') ||
    "$ENV{HOME}/.ensembl_script.conf";
  $conffile = abs_path($conffile);

  if (-e $conffile) {
    open(CONF, $conffile) or throw( 
        "Unable to open configuration file $conffile for reading: $!");

    my $serverroot = $self->serverroot;

    while (<CONF>) {
      chomp;

      # remove comments
      s/^[#;].*//;
      s/\s+[;].*$//;

      # read options into internal parameter datastructure
      next unless (/(\w\S*)\s*=\s*(.*)/);
      my $name = $1;
      my $val = $2;
      if ($val =~ /\$SERVERROOT/) {
        $val =~ s/\$SERVERROOT/$serverroot/g;
        $val = abs_path($val);
      }
      $self->param($name, $val);
    }

    $self->param('conffile', $conffile);
  }

  # override configured parameter with commandline options
  map { $self->param($_, $h{$_}) } keys %h;
  return(1);
}


=head2 parse_extra_options

  Arg[1-N]    : option descriptors that will be passed on to Getopt::Long
  Example     : $support->parse_extra_options('string_opt=s', 'numeric_opt=n');
  Description : Parse extra commandline options by passing them on to
                Getopt::Long and storing parameters in $self->param('name).
  Return type : true on success
  Exceptions  : none (caugth by $self->error)
  Caller      : general

=cut

sub parse_extra_options {
  my ($self, @params) = @_;
  
  Getopt::Long::Configure("no_pass_through");
  
  eval {
    # catch warnings to pass to $self->error
    local $SIG{__WARN__} = sub { die @_; };
    &GetOptions(\%{ $self->{'_param'} }, @params);
  };

  # comma to list
  foreach my $param (@params) {
    if ($param =~ /\@$/) {
      $param =~ s/(^\w+).*/$1/;
      $self->comma_to_list($param);
    }
  }
  
  $self->error($@) if $@;
  
  return(1);
}


=head2 allowed_params

  Arg[1-N]    : (optional) List of allowed parameters to set
  Example     : my @allowed = $self->allowed_params(qw(param1 param2));
  Description : Getter/setter for allowed parameters. This is used by
                $self->confirm_params() to avoid cluttering of output with
                conffile entries not relevant for a given script. You can use
                $self->get_common_params() as a shortcut to set them.
  Return type : Array - list of allowed parameters
  Exceptions  : none
  Caller      : general

=cut

sub allowed_params {
  my $self = shift;

  # setter
  if (@_) {
    @{ $self->{'_allowed_params'} } = @_;
  }

  # getter
  if (ref($self->{'_allowed_params'}) eq 'ARRAY') {
    return @{ $self->{'_allowed_params'} };
  } else {
    return ();
  }
}


=head2 get_common_params

  Example     : my @allowed_params = $self->get_common_params, 'extra_param';
  Description : Returns a list of commonly used parameters in the conversion
                scripts. Shortcut for setting allowed parameters with
                $self->allowed_params().
  Return type : Array - list of common parameters
  Exceptions  : none
  Caller      : general

=cut

sub get_common_params {
  return qw(
      conffile
      logpath
      logfile
      logappend
      is_component
      verbose
      interactive
      dry_run
  );
}


=head2 confirm_params

  Example     : $support->confirm_params;
  Description : Prints a table of parameters that were collected from config
                file and commandline and asks user to confirm if he wants
                to proceed.
  Return type : true on success
  Exceptions  : none
  Caller      : general

=cut

sub confirm_params {
  my $self = shift;

  if ($self->param('interactive')) {
    # print parameter table
    print "Running script with these parameters:\n\n";
    print $self->list_all_params;

    # ask user if he wants to proceed
    exit unless user_proceed("Continue?");
  }
  
  return(1);
}


=head2 list_all_params

  Example     : print LOG $support->list_all_params;
  Description : prints a table of the parameters used in the script
  Return type : String - the table to print
  Exceptions  : none
  Caller      : general

=cut

sub list_all_params {
  my $self = shift;

  $Text::Wrap::colums = 72;

  my $txt = sprintf "    %-20s%-40s\n", qw(PARAMETER VALUE);
  $txt .= "    " . "-"x70 . "\n";

  foreach my $key ($self->allowed_params) {
    my $val;
    if (defined($self->param($key))) {
      $txt .= Text::Wrap::wrap(sprintf('    %-20s', $key), ' 'x24,
        join(", ", $self->param($key)))."\n";
    }
  }

  $txt .= "\n";

  return $txt;
}


=head2 create_commandline_options

  Arg[1]      : Hashref $settings - hashref describing what to do
                Allowed keys:
                    allowed_params => 0|1   # use all allowed parameters
                    exclude => []           # listref of parameters to exclude
                    replace => {param => newval} # replace value of param with
                                                 # newval
  Example     : $support->create_commandline_options({
                    allowed_params => 1,
                    exclude => ['verbose'],
                    replace => { 'dbname' => 'homo_sapiens_vega_33_35e' }
                });
  Description : Creates a commandline options string that can be passed to any
                other script using ConversionSupport.
  Return type : String - commandline options string
  Exceptions  : none
  Caller      : general

=cut

sub create_commandline_options {
  my $self = shift;

  my ($allowed_params, $exclude, $replace) = rearrange(
    ['ALLOWED_PARAMS', 'EXCLUDE', 'REPLACE'], @_);

  my %param_hash;

  # get all allowed parameters
  if ($allowed_params) {
    # exclude params explicitly stated
    my %exclude = map { $_ => 1 } @{ $exclude || [] };

    foreach my $param ($self->allowed_params) {
      unless ($exclude{$param}) {
        my ($first, @rest) = $self->param($param);
        next unless (defined($first));

        if (@rest) {
          $first = join(",", $first, @rest);
        }
        $param_hash{$param} = $first;
      }
    }

  }

  # replace values
  foreach my $key (keys %{ $replace || {} }) {
    $param_hash{$key} = $replace->{$key};
  }

  # create the commandline options string
  my $options_string;
  foreach my $param (keys %param_hash) {
    $options_string .= sprintf("--%s %s ", $param, $param_hash{$param});
  }

  return $options_string;
}


=head2 check_required_params

  Arg[1-N]    : List @params - parameters to check
  Example     : $self->check_required_params(qw(dbname host port));
  Description : Checks $self->param to make sure the requested parameters
                have been set. Dies if parameters are missing.
  Return type : true on success
  Exceptions  : none
  Caller      : general

=cut

sub check_required_params {
  my ($self, @params) = @_;
  
  my @missing = ();
  
  foreach my $param (@params) {
    push @missing, $param unless $self->param($param);
  }
  
  if (@missing) {
    throw("Missing parameters: @missing.\nYou must specify them on the commandline or in your conffile.\n");
  }
  
  return(1);
}


=head2 comma_to_list

  Arg[1-N]    : list of parameter names to parse
  Example     : $support->comma_to_list('chromosomes');
  Description : Transparently converts comma-separated lists into arrays (to
                allow different styles of commandline options, see perldoc
                Getopt::Long for details). Parameters are converted in place
                (accessible through $self->param('name')).
  Return type : true on success
  Exceptions  : none
  Caller      : general

=cut

sub comma_to_list {
  my $self = shift;
  
  foreach my $param (@_) {
    $self->param($param, split (/,/, join (',', $self->param($param))));
  }
  
  return(1);
}


=head2 list_or_file

  Arg[1]      : Name of parameter to parse
  Example     : $support->list_or_file('gene_stable_id');
  Description : Determines whether a parameter holds a list or it is a filename
                to read the list entries from.
  Return type : true on success
  Exceptions  : thrown if list file can't be opened
  Caller      : general

=cut

sub list_or_file {
  my ($self, $param) = @_;
  
  my @vals = $self->param($param);
  return unless (@vals);

  my $firstval = $vals[0];
  
  if (scalar(@vals) == 1 && -e $firstval) {
    # we didn't get a list of values, but a file to read values from
    @vals = ();
    
    open(IN, $firstval) or throw("Cannot open $firstval for reading: $!");
    
    while(<IN>){
      chomp;
      push(@vals, $_);
    }
    
    close(IN);
    
    $self->param($param, @vals);
  }
  
  $self->comma_to_list($param);
  
  return(1);
}


=head2 param

  Arg[1]      : Parameter name
  Arg[2-N]    : (optional) List of values to set
  Example     : my $dbname = $support->param('dbname');
                $support->param('port', 3306);
                $support->param('chromosomes', 1, 6, 'X');
  Description : Getter/setter for parameters. Accepts single-value params and
                list params.
  Return type : Scalar value for single-value parameters, array of values for
                list parameters
  Exceptions  : thrown if no parameter name is supplied
  Caller      : general

=cut

sub param {
  my $self = shift;
  my $name = shift or throw("You must supply a parameter name");

  # setter
  if (@_) {
    if (scalar(@_) == 1) {
      # single value
      $self->{'_param'}->{$name} = shift;
    } else {
      # list of values
      undef $self->{'_param'}->{$name};
      @{ $self->{'_param'}->{$name} } = @_;
    }
  }

  # getter
  if (ref($self->{'_param'}->{$name}) eq 'ARRAY') {
    # list parameter
    return @{ $self->{'_param'}->{$name} };
  } elsif (defined($self->{'_param'}->{$name})) {
    # single-value parameter
    return $self->{'_param'}->{$name};
  } else {
    return undef;
  }
}


=head2 serverroot

  Arg[1]      : (optional) String - root directory of your ensembl sandbox
  Example     : my $serverroot = $support->serverroot;
  Description : Getter/setter for the root directory of your ensembl sandbox.
                This is set when ConversionSupport object is created, so
                usually only used as a getter.
  Return type : String - the server root directory
  Exceptions  : none
  Caller      : general

=cut

sub serverroot {
  my $self = shift;
  $self->{'_serverroot'} = shift if (@_);
  return $self->{'_serverroot'};
}


=head2 error

  Arg[1]      : (optional) String - error message
  Example     : $support->error("An error occurred: $@");
                exit(0) if $support->error;
  Description : Getter/setter for error messages
  Return type : String - error message
  Exceptions  : none
  Caller      : general

=cut

sub error {
    my $self = shift;
    $self->{'_error'} = shift if (@_);
    return $self->{'_error'};
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

sub logger {
  my $self = shift;

  if (@_) {
    $self->{'_logger'} = shift;

  } elsif (! $self->{'_logger'}) {
    # log to STDERR if no logger supplied
    $self->{'_logger'} = new Bio::EnsEMBL::Utils::Logger(
      -VERBOSE => $self->param('verbose'),
    );
  }

  return $self->{'_logger'};
}

1;


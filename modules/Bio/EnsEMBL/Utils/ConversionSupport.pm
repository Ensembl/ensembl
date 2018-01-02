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
  my $support = new Bio::EnsEMBL::Utils::ConversionSupport($serverroot);

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

package Bio::EnsEMBL::Utils::ConversionSupport;

use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use Text::Wrap;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use FindBin qw($Bin $Script);
use POSIX qw(strftime);
use Cwd qw(abs_path);
use DBI;
use Data::Dumper;
use Fcntl qw(:flock SEEK_END);

my $species_c = 1; #counter to be used for each database connection made

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
  (my $serverroot = shift) or throw("You must supply a serverroot.");
  my $self = {
    '_serverroot'   => $serverroot,
    '_param'        => { interactive => 1 },
    '_warnings'     => 0,
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
	       'dbname|db_name=s',
	       'host|dbhost|db_host=s',
	       'port|dbport|db_port=n',
	       'user|dbuser|db_user=s',
	       'pass|dbpass|db_pass=s',
               'prod_dbname=s',
               'prod_host=s',
               'prod_port=n',
               'prod_user=s',
               'prod_pass=s',
	       'conffile|conf=s',
	       'logfile|log=s',
               'nolog',
	       'logpath=s',
               'log_base_path=s',
	       'logappend|log_append',
	       'verbose|v',
	       'interactive|i=s',
         'hideparamlist=s',
         'joblog=s',
	       'dry_run|dry|n',
	       'help|h|?',
	     );

  # reads config file
  my $conffile = $h{'conffile'} || $self->serverroot . "/sanger-plugins/vega/conf/ini-files/Conversion.ini";
  $conffile = abs_path($conffile);
  if (-e $conffile) {
    open(my $fh, '<', $conffile) or throw( 
      "Unable to open configuration file $conffile for reading: $!");
    my $serverroot = $self->serverroot;
    while (<$fh>) {
      chomp;

      # remove comments
      s/^[#;].*//;
      s/\s+[;].*$//;

      # read options into internal parameter datastructure, removing whitespace
      next unless (/(\w\S*)\s*=\s*(\S*)\s*/);
      my $name = $1;
      my $val = $2;
      if ($val =~ /\$SERVERROOT/) {
	$val =~ s/\$SERVERROOT/$serverroot/g;
	$val = abs_path($val);
      }
      $self->param($name, $val);
    }
    close $fh;
    $self->param('conffile', $conffile);
  }
  elsif ($conffile) {
    warning("Unable to open configuration file $conffile for reading: $!");
  }

# override configured parameter with commandline options
  map { $self->param($_, $h{$_}) } keys %h;

  return (1) if $self->param('nolog');

  # if logpath & logfile are not set, set them here to /ensemblweb/vega_dev/shared/logs/conversion/DBNAME/SCRIPNAME_NN.log
  if (! defined($self->param('log_base_path')))  {
    $self->param('log_base_path','/ensemblweb/shared/logs/conversion/');
  }
  my $dbname = $self->param('dbname');
  $dbname =~ s/^vega_//;
  if (not (defined($self->param('logpath')) )){
    $self->param('logpath', $self->param('log_base_path')."/".$dbname."/" );
  }
  if ( not defined $self->param('logfile') ){
    my $log = $Script;
    $log =~ s/.pl$//g;
    my $counter;
    for ($counter=1 ; (-e $self->param('logpath')."/".$log."_".sprintf("%03d", $counter).".log"); $counter++){
      #      warn  $self->param('logpath')."/".$log."_".$counter.".log";
    }
    $self->param('logfile', $log."_".sprintf("%03d", $counter).".log");
  }
  return(1);
}

=head2 parse_extra_options

  Arg[1-N]    : option descriptors that will be passed on to Getopt::Long
  Example     : $support->parse_extra_options('string_opt=s', 'numeric_opt=n');
  Description : Parse extra commandline options by passing them on to
                Getopt::Long and storing parameters in $self->param('name).
                If the last option is '...' then unknown options are kept
                in @ARGV for future calls to this method. Without '...'
                then it is an error to have further options.
  Return type : true on success
  Exceptions  : none (caugth by $self->error)
  Caller      : general

=cut

sub parse_extra_options {
  my ($self, @params) = @_;
  if(@params and $params[-1] eq "...") {
    pop @params;
    Getopt::Long::Configure("pass_through");
  } else {
    Getopt::Long::Configure("no_pass_through");
  }
  eval {
    # catch warnings to pass to $self->error
    local $SIG{__WARN__} = sub { die @_; };
    &GetOptions(\%{ $self->{'_param'} }, @params);
  };
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
	    dbname
	    host
	    port
	    user
	    pass
            nolog
	    logpath
            log_base_path
	    logfile
	    logappend
	    verbose
	    interactive
      hideparamlist
      joblog
	    dry_run
	  );
}

=head2 get_loutre_params

  Arg         : (optional) return a list to parse or not
  Example     : $support->parse_extra_options($support->get_loutre_params('parse'))
  Description : Returns a list of commonly used loutre db parameters - parse option is
                simply used to distinguish between reporting and parsing parameters
  Return type : Array - list of common parameters
  Exceptions  : none
  Caller      : general

=cut

sub get_loutre_params {
  my ($self,$p) = @_;
  if ($p) {
    return qw(
	      loutrehost=s
	      loutreport=s
	      loutreuser=s
	      loutrepass:s
	      loutredbname=s
	    );
  }
  else {
    return qw(
	      loutrehost
	      loutreport
	      loutreuser
	      loutrepass
	      loutredbname
	    );
  }
}

=head2 get_annotrack_params

  Arg         : (optional) return a list to parse or not
  Example     : $support->parse_extra_options($support->get_annotrack_params('parse'))
  Description : Returns a list of parameters used to connect to annotrack
  Return type : Array - list of common parameters
  Exceptions  : none
  Caller      : general

=cut

sub get_annotrack_params {
  my ($self,$p) = @_;
  if ($p) {
    return qw(
              annotrackhost=s
	      annotrackport=s
	      annotrackuser=s
	      annotrackpass=s
	      annotrackdbname=s
              ignoreannotrack
	    );
  }
  else {
    return qw(
	      annotrackhost
	      annotrackport
	      annotrackuser
	      annotrackpass
	      annotrackdbname
              ignoreannotrack
	    );
  }
}

=head2 remove_vega_params

  Example     : $support->remove_vega_params;
  Description : Removes Vega db conection parameters. Usefull to avoid clutter in log files when
                working exclusively with loutre
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub remove_vega_params {
  my $self = shift;
  foreach my $param (qw(dbname host port user pass)) {
    $self->{'_param'}{$param} = undef;
  }
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

  # print parameter table
  return 1 if($self->param('hideparamlist') and not $self->param('interactive'));
  print "Running script with these parameters:\n\n";
  print $self->list_all_params;

  if ($self->param('host') eq 'ensweb-1-10') {
    # ask user if he wants to proceed
    exit unless $self->user_proceed("**************\n\n You're working on ensdb-1-10! Is that correct and you want to continue ?\n\n**************");
  }
  else {
    # ask user if he wants to proceed
    exit unless $self->user_proceed("Continue?");
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
  my $txt = sprintf "    %-25s%-110s\n", qw(PARAMETER VALUE);
  $txt .= "    " . "-"x136 . "\n";
  $Text::Wrap::columns = 130;
  my @params = $self->allowed_params;
  foreach my $key (@params) {
    my @vals = $self->param($key);
    if (@vals) {
      $txt .= Text::Wrap::wrap( sprintf('   %-25s', $key),
				' 'x28,
				join(", ", @vals)
			      ) . "\n";
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
  my ($self, $settings, $param_hash) = @_;
  my %param_hash = $param_hash ? %$param_hash : ();

  # get all allowed parameters
  if ($settings->{'allowed_params'}) {
    # exclude params explicitly stated
    my %exclude = map { $_ => 1 } @{ $settings->{'exclude'} || [] };
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
  foreach my $key (keys %{ $settings->{'replace'} || {} }) {
    $param_hash{$key} = $settings->{'replace'}->{$key};
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
    push @missing, $param unless defined $self->param($param);
  }
  if (@missing) {
    throw("Missing parameters: @missing.\nYou must specify them on the commandline or in your conffile.\n");
  }
  return(1);
}

=head2 user_proceed

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub user_proceed {
  my ($self, $text) = @_;

  if ($self->param('interactive')) {
    print "$text\n" if $text;
    print "[y/N] ";
    my $input = lc(<>);
    chomp $input;
    unless ($input eq 'y') {
      print "Skipping.\n";
      return(0);
    }
  }

  return(1);
}


=head2 read_user_input

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : my $ret = $support->read_user_input("Choose a number [1/2/3]");
                if ($ret == 1) {
                    # do something
                } elsif ($ret == 2) {
                    # do something else
                }
  Description : If running interactively, the user is asked for input.
  Return type : String - user's input
  Exceptions  : none
  Caller      : general

=cut

sub read_user_input {
  my ($self, $text) = @_;

  if ($self->param('interactive')) {
    print "$text\n" if $text;
    my $input = <>;
    chomp $input;
    return $input;
  }
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
    $self->param($param,
		 split (/,/, join (',', $self->param($param))));
  }
  return(1);
}

=head2 list_or_file

  Arg[1]      : Name of parameter to parse
  Example     : $support->list_or_file('gene');
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
    open(my $fh, '<', $firstval) or throw("Cannot open $firstval for reading: $!");
    while(<$fh>){
      chomp;
      push(@vals, $_);
    }
    close($fh);
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
    return ();
  }
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

=head2 get_database

  Arg[1]      : String $database - the type of database to connect to
                (eg core, otter)
  Arg[2]      : (optional) String $prefix - the prefix used for retrieving the
                connection settings from the configuration
  Example     : my $db = $support->get_database('core');
  Description : Connects to the database specified.
  Return type : DBAdaptor of the appropriate type
  Exceptions  : thrown if asking for unknown database
  Caller      : general

=cut

sub get_database {
  my $self = shift;
  my $database = shift or throw("You must provide a database");
  my $prefix = shift || '';
  $self->check_required_params(
    "${prefix}host",
    "${prefix}port",
    "${prefix}user",
    "${prefix}dbname",
  );
  my %adaptors = (
    core    => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    ensembl => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    evega   => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    otter   => 'Bio::Otter::DBSQL::DBAdaptor',
    vega    => 'Bio::Otter::DBSQL::DBAdaptor',
    compara => 'Bio::EnsEMBL::Compara::DBSQL::DBAdaptor',
    loutre  => 'Bio::Vega::DBSQL::DBAdaptor',
  );
  throw("Unknown database: $database") unless $adaptors{$database};

  $self->dynamic_use($adaptors{$database});
  my $species = 'species' . $species_c;
  my $dba = $adaptors{$database}->new(
    -host    => $self->param("${prefix}host"),
    -port    => $self->param("${prefix}port"),
    -user    => $self->param("${prefix}user"),
    -pass    => $self->param("${prefix}pass") || '',
    -dbname  => $self->param("${prefix}dbname"),
    -group   => $database,
    -species => $species,
  );
  #can use this approach to get dna from another db
#  my $dna_db = $adaptors{$database}->new(
#    -host => 'otterlive',
#    -port => '3301',
#    -user => $self->param("${prefix}user"),
#    -pass => $self->param("${prefix}pass"),
#    -dbname => 'loutre_human',
#  );
#  $dba->dnadb($dna_db);

  # otherwise explicitely set the dnadb to itself - by default the Registry assumes
  # a group 'core' for this now
  $dba->dnadb($dba);

  $species_c++;

  $self->{'_dba'}->{$database} = $dba;
  $self->{'_dba'}->{'default'} = $dba unless $self->{'_dba'}->{'default'};
  return $self->{'_dba'}->{$database};
}


=head2 get_dbconnection

  Arg[1]      : (optional) String $prefix - the prefix used for retrieving the
                connection settings from the configuration
  Example     : my $dbh = $self->get_dbconnection;
  Description : Connects to the database server specified. You don't have to
                specify a database name (this is useful for running commands
                like $dbh->do('show databases')).
  Return type : DBI database handle
  Exceptions  : thrown if connection fails
  Caller      : general
  Status      : At Risk

=cut

sub get_dbconnection {
  my $self = shift;
  my $prefix = shift;
 
  $self->check_required_params(
      "${prefix}host",
      "${prefix}port",
      "${prefix}user",
  );

  my $dsn = "DBI:" . ($self->param('driver')||'mysql') .
            ":host=" . $self->param("${prefix}host") .
            ";port=" . $self->param("${prefix}port");

  if ($self->param("${prefix}dbname")) {
    $dsn .= ";dbname=".$self->param("${prefix}dbname");
  }

#  warn $dsn;

  my $dbh;
  my $user = $self->param("${prefix}user");
  my $pass = $self->param("${prefix}pass");
  eval{
    $dbh = DBI->connect($dsn, $user, $pass, {'RaiseError' => 1, 'PrintError' => 0});
  };

  if (!$dbh || $@ || !$dbh->ping) {
    $self->log_error("Could not connect to db server as user ".
      $self->param("${prefix}user") .
      " using [$dsn] as a locator:\n" . $DBI::errstr . $@);
  }

  $self->{'_dbh'} = $dbh;
  return $self->{'_dbh'};

}


=head2 dba

  Arg[1]      : (optional) String $database - type of db apaptor to retrieve
  Example     : my $dba = $support->dba;
  Description : Getter for database adaptor. Returns default (i.e. created
                first) db adaptor if no argument is provided.
  Return type : Bio::EnsEMBL::DBSQL::DBAdaptor or Bio::Otter::DBSQL::DBAdaptor
  Exceptions  : none
  Caller      : general

=cut

sub dba {
  my ($self, $database) = @_;
  return $self->{'_dba'}->{$database} || $self->{'_dba'}->{'default'};
}

=head2 dynamic_use

  Arg [1]    : String $classname - The name of the class to require/import
  Example    : $self->dynamic_use('Bio::EnsEMBL::DBSQL::DBAdaptor');
  Description: Requires and imports the methods for the classname provided,
               checks the symbol table so that it doesnot re-require modules
               that have already been required.
  Returntype : true on success
  Exceptions : Warns to standard error if module fails to compile
  Caller     : internal

=cut

sub dynamic_use {
  my ($self, $classname) = @_;
  my ($parent_namespace, $module) = $classname =~/^(.*::)(.*)$/ ? ($1,$2) : ('::', $classname);

  no strict 'refs'; ## no critic
  # return if module has already been imported
  return 1 if $parent_namespace->{$module.'::'} && %{ $parent_namespace->{$module.'::'}||{} };

  eval "require $classname"; ## no critic
  throw("Failed to require $classname: $@") if ($@);
  $classname->import();

  return 1;
}

=head2 get_chrlength

  Arg[1]      : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor $dba
  Arg[2]      : (optional) String $version - coord_system version
  Arg[3]      : (optional) String $type - type of region eg chromsome (defaults to 'toplevel')
  Arg[4]      : (optional) Boolean - return non reference slies as well (required for haplotypes eq 6-COX)
  Arg[5]      : (optional) Override chromosome parameter filtering with this array reference. Empty denotes all.
  Example     : my $chr_length = $support->get_chrlength($dba);
  Description : Get all chromosomes and their length from the database. Return
                chr_name/length for the chromosomes the user requested (or all
                chromosomes by default)
  Return type : Hashref - chromosome_name => length
  Exceptions  : thrown if not passing a Bio::EnsEMBL::DBSQL::DBAdaptor
  Caller      : general

=cut

sub get_chrlength {
  my ($self, $dba, $version,$type,$include_non_reference,$chroms) = @_;
  $dba  ||= $self->dba;
  $type ||= 'toplevel';

  throw("get_chrlength should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n")
    unless ($dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));

  my $sa = $dba->get_SliceAdaptor;

  my @chromosomes = map { $_->seq_region_name } 
    @{ $sa->fetch_all($type,$version,$include_non_reference) };
  my %chr = map { $_ => $sa->fetch_by_region($type, $_, undef, undef, undef, $version)->length } @chromosomes;

  my @wanted = $self->param('chromosomes');
  @wanted = @$chroms if defined $chroms and ref($chroms) eq 'ARRAY';

  if (@wanted) {
    # check if user supplied invalid chromosome names
    foreach my $chr (@wanted) {
      my $found = 0;
      foreach my $chr_from_db (keys %chr) {
	if ($chr_from_db eq $chr) {
	  $found = 1;
	  last;
	}
      }
      unless ($found) {
        #if you get this when using ensembl-vega then the dbname will be misleading.
	warning("Didn't find chromosome $chr in database " .
		  $self->param('dbname'));
      }
    }

    # filter to requested chromosomes only
  HASH:
    foreach my $chr_from_db (keys %chr) {
      foreach my $chr (@wanted) {
	if ($chr_from_db eq $chr) {
	  next HASH;
	}
      }
      delete($chr{$chr_from_db});
    }
  }

  return \%chr;
}

=head2 get_ensembl_chr_mapping

  Arg[1]      : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor $dba
  Arg[2]      : (optional) String $version - coord_system version
  Example     : my $ensembl_mapping = $support->get_ensembl_chr_mapping($dba);
  Description : Gets a mapping between Vega chromosome names and their
                equivalent Ensembl chromosomes. Works with non-reference chromosomes
  Return type : Hashref - Vega name => Ensembl name
  Exceptions  : thrown if not passing a Bio::EnsEMBL::DBSQL::DBAdaptor
  Caller      : general

=cut

sub get_ensembl_chr_mapping {
  my ($self, $dba, $version) = @_;
  $dba ||= $self->dba;
  throw("get_ensembl_chr_mapping should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n") unless ($dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));

  my $sa = $dba->get_SliceAdaptor;
  my @chromosomes = map { $_->seq_region_name } 
    @{ $sa->fetch_all('chromosome', $version, 1) };

  my %chrs;
  foreach my $chr (@chromosomes) {
    my $sr = $sa->fetch_by_region('chromosome', $chr, undef, undef, undef, $version);
    my ($ensembl_name_attr) = @{ $sr->get_all_Attributes('ensembl_name') };
    if ($ensembl_name_attr) {
      $chrs{$chr} = $ensembl_name_attr->value;
    } else {
      $chrs{$chr} = $chr;
    }
  }
  return \%chrs;
}

=head2 get_taxonomy_id

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBAdaptor $dba
  Example     : my $sid = $support->get_taxonony_id($dba);
  Description : Retrieves the taxononmy ID from the meta table
  Return type : Int - the taxonomy ID
  Exceptions  : thrown if no taxonomy ID is found in the database
  Caller      : general

=cut

sub get_taxonomy_id {
  my ($self, $dba) = @_;
  $dba ||= $self->dba;
  my $sql = 'SELECT meta_value FROM meta WHERE meta_key = "species.taxonomy_id"';
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute;
  my ($tid) = $sth->fetchrow_array;
  $sth->finish;
  $self->throw("Could not determine taxonomy_id from database.") unless $tid;
  return $tid;
}

=head2 get_species_scientific_name

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBAdaptor $dba
  Example     : my $species = $support->get_species_scientific_name($dba);
  Description : Retrieves the species scientific name (Genus species) from the
                meta table
  Return type : String - species scientific name
  Exceptions  : thrown if species name can not be determined from db
  Caller      : general

=cut

sub get_species_scientific_name {
  my ($self, $dba) = @_;
  $dba ||= $self->dba;
  my $sql = "SELECT meta_value FROM meta WHERE meta_key = \'species.scientific_name\'";
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute;
  my @sp;
  while (my @row = $sth->fetchrow_array) {
    push @sp, $row[0];
  }
  if (! @sp || @sp > 1) {
    $self->throw("Could not retrieve a single species scientific name from database.");
  }
  $sth->finish;
  my $species = $sp[0];
  $species =~ s/ /_/g;
  return $species;
}

=head2 species

  Arg[1]      : (optional) String $species - species name to set
  Example     : my $species = $support->species;
                my $url = "http://vega.sanger.ac.uk/$species/";
  Description : Getter/setter for species name (Genus_species). If not set, it's
                determined from database's meta table
  Return type : String - species name
  Exceptions  : none
  Caller      : general

=cut

sub species {
  my $self = shift;
  $self->{'_species'} = shift if (@_);
  # get species name from database if not set
  unless ($self->{'_species'}) {
    $self->{'_species'} = $self->get_species_scientific_name;
  }
  return $self->{'_species'};
}

=head2 sort_chromosomes

  Arg[1]      : (optional) Hashref $chr_hashref - Hashref with chr_name as keys
  Example     : my $chr = { '6-COX' => 1, '1' => 1, 'X' => 1 };
                my @sorted = $support->sort_chromosomes($chr);
  Description : Sorts chromosomes in an intuitive way (numerically, then
                alphabetically). If no chromosome hashref is passed, it's
                retrieve by calling $self->get_chrlength()
  Return type : List - sorted chromosome names
  Exceptions  : thrown if no hashref is provided
  Caller      : general

=cut

sub sort_chromosomes {
  my ($self, $chr_hashref, $version, $include_non_reference) = @_;
  $chr_hashref = $self->get_chrlength($self->dba, $version, 'chromosome', $include_non_reference) unless ($chr_hashref);
  throw("You have to pass a hashref of your chromosomes")
    unless ($chr_hashref and ref($chr_hashref) eq 'HASH');
  return (sort _by_chr_num keys %$chr_hashref);
}

=head2 _by_chr_num

  Example     : my @sorted = sort _by_chr_num qw(X, 6-COX, 14, 7);
  Description : Subroutine to use in sort for sorting chromosomes. Sorts
                numerically, then alphabetically
  Return type : values to be used by sort
  Exceptions  : none
  Caller      : internal ($self->sort_chromosomes)

=cut

sub _by_chr_num {
  my @awords = split /-/, $a;
  my @bwords = split /-/, $b;

  my $anum = $awords[0];
  my $bnum = $bwords[0];

  if ($anum !~ /^[0-9]*$/) {
    if ($bnum !~ /^[0-9]*$/) {
      return $anum cmp $bnum;
    } else {
      return 1;
    }
  }
  if ($bnum !~ /^[0-9]*$/) {
    return -1;
  }

  if ($anum <=> $bnum) {
    return $anum <=> $bnum;
  } else {
    if ($#awords == 0) {
      return -1;
    } elsif ($#bwords == 0) {
      return 1;
    } else {
      return $awords[1] cmp $bwords[1];
    }
  }
}

=head2 split_chromosomes_by_size

  Arg[1]      : (optional) Int $cutoff - the cutoff in bp between small and
                large chromosomes
  Arg[2]      : (optional) Boolean to include duplicate regions, ie PAR or not
                (default is no)
  Arg[3]      : (optional) Coordsystem version to retrieve

  Example     : my $chr_slices = $support->split_chromosomes_by_size;
                foreach my $block_size (keys %{ $chr_slices }) {
                    print "Chromosomes with blocksize $block_size: ";
                    print join(", ", map { $_->seq_region_name }
                                        @{ $chr_slices->{$block_size} });
                }
  Description : Determines block sizes for storing DensityFeatures on
                chromosomes, and return slices for each chromosome. The block
                size is determined so that you have 150 bins for the smallest
                chromosome over 5 Mb in length. For chromosomes smaller than 5
                Mb, an additional smaller block size is used to yield 150 bins
                for the overall smallest chromosome. This will result in
                reasonable resolution for small chromosomes and high
                performance for big ones. Does not return non-reference seq_regions
  Return type : Hashref (key: block size; value: Arrayref of chromosome
                Bio::EnsEMBL::Slices)
  Exceptions  : none
  Caller      : density scripts

=cut

sub split_chromosomes_by_size {
  my $self   = shift;
  my $cutoff = shift || 5000000;
  my $dup    = shift || 0;
  my $cs_version = shift;
  my $include_non_reference = 1; #get non reference slices
  my $slice_adaptor = $self->dba->get_SliceAdaptor;
  my ($top_slices, $wanted_slices, $missed_slices);
  if ($self->param('chromosomes')) {
    foreach my $chr ($self->param('chromosomes')) {
      push @{ $top_slices }, $slice_adaptor->fetch_by_region('chromosome', $chr);
    }
  }
  else {
    $top_slices = $slice_adaptor->fetch_all('chromosome',$cs_version,$include_non_reference,$dup);
  }

  # filter out patches, if present
  $wanted_slices = [ grep { $_->is_reference or $self->is_haplotype($_,$self->dba) } @$top_slices ];

  # make a note of which ones are excluded
  $missed_slices = [ grep { ! $_->is_reference and ! $self->is_haplotype($_,$self->dba) } @$top_slices ];

  #warn non PATCH toplevels slice that are excluded (could be haplotypes if an earlier stage has failed)
  foreach (@{$missed_slices || []}) {
    my $sr_name = $_->seq_region_name;
    if ($self->is_patch($_)) {
      $self->log("Excluding $sr_name\n");
    }
    else {
      $self->log_warning("Excluding $sr_name, are you sure ?\n");
    }
  }

  my ($big_chr, $small_chr, $min_big_chr, $min_small_chr);
  foreach my $slice (@{ $wanted_slices }) {
    next if ($slice->length eq 10000); #hack for chrY pseudoslice
    if ($slice->length < $cutoff) {
      if (! $min_small_chr or ($min_small_chr > $slice->length)) {
	$min_small_chr = $slice->length;
      }
      # push small chromosomes onto $small_chr
      push @{ $small_chr }, $slice;
    }
    elsif (! $min_big_chr or ($min_big_chr > $slice->length) ){
      $min_big_chr = $slice->length;
    }
    # push _all_ chromosomes onto $big_chr
    push @{ $big_chr }, $slice;
  }
  my $chr_slices;
  $chr_slices->{int($min_big_chr/150)}   = $big_chr if $min_big_chr;
  $chr_slices->{int($min_small_chr/150)} = $small_chr if $min_small_chr;
  return $chr_slices;
}

=head2 is_patch

  Arg[1]      : B::E::Slice
  Example     : if ($support->is_patch($slice)) { ...
  Description : Looks at seq_region attributes to decide if a slice is a patch or not
                If PATCH seq_region_attrib is not there check to see if name suggests it is a PATCH
  Return type : true/false

=cut

sub is_patch {
  my ($self,$slice) = @_;
  my @patch_attrib_types = qw(patch_fix patch_novel); #seq_region_attribs used to define a patch
  foreach my $attrib_type (@patch_attrib_types) {
    if (@{$slice->get_all_Attributes($attrib_type)}) {
      return 1;
    }
  }
  return 0;
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

  # strip off leading linebreaks so that indenting doesn't break
  $txt =~ s/^(\n*)//;

  $txt = $1."    "x$indent . $txt;
  my $fh = $self->{'_log_filehandle'};
  throw("Unable to obtain log filehandle") unless $fh;
  print $fh "$txt";
  return(1);
}

=head2 lock_log

  Description : Use flock-style locks to lock log and fastforward to end.
                Useful if log is being written to by multiple processes.
=cut

sub lock_log {
  my ($self) = @_;
  ## no critic
  my $fh = $self->{'_log_filehandle'};
  return if -t $fh or -p $fh; # Shouldn't lock such things
  flock($self->{'_log_filehandle'},LOCK_EX) || return 0;
  seek($self->{'_log_filehandle'},0,SEEK_END); # fail ok, prob not reg file
  return 1;
}

=head2 unlock_log

  Description : Unlock log previously locked by lock_log.

=cut

sub unlock_log {
  my ($self) = @_;
  ## no critic
  my $fh = $self->{'_log_filehandle'};
  return if -t $fh or -p $fh; # We don't lock such things
  # flush is implicit in flock
  flock($self->{'_log_filehandle'},LOCK_UN) || return 0;
  return 1;
}

=head2 log_warning

  Arg[1]      : String $txt - the warning text to log
  Arg[2]      : Int $indent - indentation level for log message
  Arg[3]      : Bool - add a line break before warning if true
  Example     : my $log = $support->log_filehandle;
                $support->log_warning('Log foo.\n', 1);
  Description : Logs a message via $self->log and increases the warning counter.
  Return type : true on success
  Exceptions  : none
  Caller      : general

=cut

sub log_warning {
  my ($self, $txt, $indent, $break) = @_;
  $txt = "WARNING: " . $txt;
  $txt = "\n$txt" if ($break);
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
  $self->log("Exiting.\n");
  exit;
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
  return(0) unless $self->param('verbose');
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
  my ($self, $mode, $date) = @_;
  $mode ||= '>';
  $mode = '>>' if ($self->param('logappend'));
  my $fh = \*STDERR;
  if (my $logfile = $self->param('logfile')) {
    $logfile .= "_$date" if $date;
    if (my $logpath = $self->param('logpath')) {
      unless (-e $logpath) {
	system("mkdir $logpath") == 0 or
	  $self->log_error("Can't create log dir $logpath: $!\n");
      }
      $logfile = "$logpath/$logfile";
    }
    open($fh, "$mode", $logfile) or throw(
      "Unable to open $logfile for writing: $!");
  }
  $self->{'_log_filehandle'} = $fh;
  return $self->{'_log_filehandle'};
}

=head2 filehandle

  Arg[1]      : String $mode - file access mode
  Arg[2]      : String $file - input or output file
  Example     : my $fh = $support->filehandle('>>', '/path/to/file');
                # print to the filehandle
                print $fh 'Your text goes here...\n';
  Description : Returns a filehandle (*STDOUT for writing, *STDIN for reading
                by default) to print to or read from.
  Return type : Filehandle - the filehandle
  Exceptions  : thrown if file can't be opened
  Caller      : general

=cut

sub filehandle {
  my ($self, $mode, $file) = @_;
  $mode ||= ">";
  my $fh;
  if ($file) {
    open($fh, "$mode", $file) or throw(
      "Unable to open $file for writing: $!");
  } elsif ($mode =~ />/) {
    $fh = \*STDOUT;
  } elsif ($mode =~ /</) {
    $fh = \*STDIN;
  }
  return $fh;
}

=head2 init_log_date

  Example     : $support->init_log_date;
  Description : Opens a filehandle to a logfile with the date in the file name
  Return type : Filehandle - the log filehandle
  Exceptions  : none
  Caller      : general

=cut

sub init_log_date {
  my $self = shift;
  my $date = $self->date;
  return $self->init_log($date);
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
  my $date = shift;

  # get a log filehandle
  my $log = $self->log_filehandle(undef,$date);

  # print script name, date, user who is running it
  unless($self->param('hideparamlist')) {
    my $hostname = `hostname`;
    chomp $hostname;
    my $script = "$hostname:$Bin/$Script";
    my $user = `whoami`;
    chomp $user;
    $self->log("Script: $script\nDate: ".$self->date_and_time."\nUser: $user\n");

    # print parameters the script is running with
    $self->log("Parameters:\n\n");
    $self->log($self->list_all_params);
  }
  # remember start time
  $self->{'_start_time'} = time;

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

  unless($self->param('hideparamlist')) {
    $self->log("\nAll done. ".$self->warnings." warnings. ");
    if ($self->{'_start_time'}) {
      $self->log("Runtime ");
      my $diff = time - $self->{'_start_time'};
      my $sec = $diff % 60;
      $diff = ($diff - $sec) / 60;
      my $min = $diff % 60;
      my $hours = ($diff - $min) / 60;
      $self->log("${hours}h ${min}min ${sec}sec ");
    }
    $self->log($self->date_and_mem."\n\n");
  }
  if($self->param('joblog')) {
    my $fh;
    unless(open($fh,'>',$self->param('joblog'))) {
      $self->log_warning("Could not log job to '".$self->param('joblog').
                         "': $!");
      return 1;
    }
    my @keys;
    push @keys,"RUNTIME",(time - $self->{'_start_time'});
    my $mem = `ps -p $$ -o vsz |tail -1`;
    chomp $mem;
    push @keys,"MEMORY",$mem;
    push @keys,"WARNINGS",$self->warnings;
    push @keys,"START",$self->{'_start_time'};
    push @keys,"END",time;
    while(@keys) {
      my ($k,$v) = splice(@keys,0,2);
      print $fh "$k: $v\n"; 
    }
    close $fh;
  }
  return(1);
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
  Description : Prints a nicely formatted datetamp (YYYY-DD-MM)
  Return type : String - the timestamp
  Exceptions  : none
  Caller      : general

=cut

sub date {
  return strftime "%Y-%m-%d", localtime;
}

=head2 date_and_time

  Example     : print "Date: " . $support->date . "\n";
  Description : Prints a nicely formatted timestamp (YYYY-DD-MM hh:mm:ss)
  Return type : String - the timestamp
  Exceptions  : none
  Caller      : general

=cut

sub date_and_time {
  return strftime "%Y-%m-%d %T", localtime;
}

=head2 format_time

  Example     : print $support->format_time($gene->modifed_date) . "\n";
  Description : Prints timestamps from the database
  Return type : String - nicely formatted time stamp
  Exceptions  : none
  Caller      : general

=cut


sub date_format {
  my( $self, $time, $format ) = @_;
  my( $d,$m,$y) = (localtime($time))[3,4,5];
  my %S = ('d'=>sprintf('%02d',$d),'m'=>sprintf('%02d',$m+1),'y'=>$y+1900);
  (my $res = $format ) =~s/%(\w)/$S{$1}/ge;
  return $res;
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

=head2 commify

  Arg[1]      : Int $num - a number to commify
  Example     : print "An easy to read number: ".$self->commify(100000000);
                # will print 100,000,000
  Description : put commas into a number to make it easier to read
  Return type : a string representing the commified number
  Exceptions  : none
  Caller      : general
  Status      : stable

=cut

sub commify {
  my $self = shift;
  my $num = shift;

  $num = reverse($num);
  $num =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;

  return scalar reverse $num;
}

=head2 fetch_non_hidden_slices

  Arg[1]      : B::E::SliceAdaptor
  Arg[2]      : B::E::AttributeAdaptor
  Arg[3]      : string $coord_system_name (optional) - 'chromosome' by default
  Arg[4]      : string $coord_system_version (optional) - 'otter' by default
  Example     : $chroms = $support->fetch_non_hidden_slice($sa,$aa);
  Description : retrieve all slices from a loutre database that don't have a hidden attribute.
                Doesn't retrieve non-reference slices
  Return type : arrayref
  Caller      : general
  Status      : stable

=cut

sub fetch_non_hidden_slices {
  my $self = shift;
  my $aa   = shift or throw("You must supply an attribute adaptor");
  my $sa   = shift or throw("You must supply a slice adaptor");
  my $cs   = shift || 'chromosome';
  my $cv   = shift || 'Otter';
  my $visible_chroms;
  foreach my $chrom ( @{$sa->fetch_all($cs,$cv)} ) {
    my $chrom_name = $chrom->name;
    my $attribs = $aa->fetch_all_by_Slice($chrom,'hidden');
    if ( scalar(@$attribs) > 1 ) {
      $self->log_warning("More than one hidden attribute for chromosome $chrom_name\n");
    }
    elsif ($attribs->[0]->value == 0) {				
      push @$visible_chroms, $chrom;
    }
    elsif ($attribs->[0]->value == 1) {	
      $self->log_verbose("chromosome $chrom_name is hidden\n");	
    }
    else {
      $self->log_warning("No hidden attribute for chromosome $chrom_name\n");
    }
  }
  return $visible_chroms;
}

=head2 get_non_hidden_slice_names

  Arg[1]      : B::E::SliceAdaptor
  Arg[2]      : B::E::AttributeAdaptor
  Arg[3]      : string $coord_system_name (optional) - 'chromosome' by default
  Arg[4]      : string $coord_system_version (optional) - 'otter' by default
  Example     : $chrom_names = $support->get_non_hidden_slice_names($sa,$aa);
  Description : retrieve names of all slices from a loutre database that don't have a hidden attribute.
                Doesn't retrieve non-reference slices
  Return type : arrayref of names of all non-hidden slices
  Caller      : general
  Status      : stable

=cut

sub get_non_hidden_slice_names {
  my $self = shift;
  my $aa   = shift or throw("You must supply an attribute adaptor");
  my $sa   = shift or throw("You must supply a slice adaptor");
  my $cs   = shift || 'chromosome';
  my $cv   = shift || 'Otter';
  my $visible_chrom_names;
  foreach my $chrom ( @{$sa->fetch_all($cs,$cv)} ) {
    my $chrom_name = $chrom->seq_region_name;
    my $attribs = $aa->fetch_all_by_Slice($chrom,'hidden');
    if ( scalar(@$attribs) > 1 ) {
      $self->log_warning("More than one hidden attribute for chromosome $chrom_name\n");
    }
    elsif ($attribs->[0]->value == 0) {				
      push @$visible_chrom_names, $chrom_name;
    }
    elsif ($attribs->[0]->value == 1) {	
      $self->log_verbose("chromosome $chrom_name is hidden\n");	
    }
    else {
      $self->log_warning("No hidden attribute for chromosome $chrom_name\n");
    }
  }
  return $visible_chrom_names;
}


=head2 get_wanted_chromosomes

  Arg[1]      : B::E::SliceAdaptor
  Arg[2]      : B::E::AttributeAdaptor
  Arg[3]      : string $coord_system_name (optional) - 'chromosome' by default
  Arg[4]      : string $coord_system_version (optional) - 'otter' by default
  Example     : $chr_names = $support->get_wanted_chromosomes($laa,$lsa);
  Description : retrieve names of slices from a lutra database that are ready for dumping to Vega.
                Deals with list of names to ignore (ignore_chr = LIST)
  Return type : arrayref of slices
  Caller      : general
  Status      : stable

=cut

sub get_wanted_chromosomes {
  my $self = shift;
  my $aa   = shift or throw("You must supply an attribute adaptor");
  my $sa   = shift or throw("You must supply a slice adaptor");
  my $cs   = shift || 'chromosome';
  my $cv   = shift || 'Otter';
  my $export_mode = $self->param('release_type');
  my $release = $self->param('vega_release');
  my $names;
  my $chroms  = $self->fetch_non_hidden_slices($aa,$sa,$cs,$cv); 
  CHROM:
  foreach my $chrom (@$chroms) {
    my $attribs = $aa->fetch_all_by_Slice($chrom);
    my $vals = $self->get_attrib_values($attribs,'vega_export_mod');
    if (scalar(@$vals > 1)) {
      $self->log_warning ("Multiple attribs for \'vega_export_mod\', please fix before continuing");
      exit;
    }
    next CHROM if (! grep { $_ eq $export_mode} @$vals);
    $vals =  $self->get_attrib_values($attribs,'vega_release',$release);	
    if (scalar(@$vals > 1)) {
      $self->log_warning ("Multiple attribs for \'vega_release\' value = $release , please fix before continuing");
      exit;
    }
    next CHROM if (! grep { $_ eq $release} @$vals);
    my $name = $chrom->seq_region_name;
    if (my @ignored = $self->param('ignore_chr')) {
      next CHROM if (grep {$_ eq $name} @ignored);
    }
    push @{$names}, $name;
  }
  return $names;
}

=head2 is_haplotype

  Arg[1]      : B::E::Slice
  Arg[2]:     : B::E::DBAdaptor (optional, if you don't supply one then the *first* one you generated is returned, which may or may not be what you want!)
  Description : Is the slice a Vega haplotype? At the moment this is 
    implemented by testing for presence of vega_ref_chrom but non_ref
    which is correct in practice, but really misses the prupose of
    vega_ref_chrom, so this might bite us if that changes.
  Return type : boolean

=cut

sub is_haplotype {
  my ($self,$slice,$dba) = @_;

  $dba ||= $self->dba;
  my $aa = $dba->get_adaptor('Attribute');

  my $attribs = $aa->fetch_all_by_Slice($slice);
  return (@{$self->get_attrib_values($attribs,'vega_ref_chrom')} and
          @{$self->get_attrib_values($attribs,'non_ref',1)});
}

=head2 get_unique_genes

  Arg[1]      : B::E::Slice
  Arg[2]      : B::E::DBAdaptor (optional, if you don't supply one then the *first* one you generated is returned, which may or may not be what you want!)
  Example     : $genes = $support->get_unique_genes($slice,$dba);
  Description : Retrieve genes that are only on the slice itself - used for human where assembly patches
                are in the assembly_exception table. Needs the PATCHes to have 'non_ref' seq_region_attributes.
  Return type : arrayref of genes
  Caller      : general
  Status      : stable

=cut

sub get_unique_genes {
  my $self  = shift;
  my ($slice,$dba) = @_;
  $slice or throw("You must supply a slice");
  $dba ||= $self->dba;

  my $sa    = $dba->get_adaptor('Slice');
  my $ga    = $dba->get_adaptor('Gene');
  my $patch = 0;
  my $genes = [];
  if ( ! $slice->is_reference() and ! $self->is_haplotype($slice,$dba) ) {
#  if ( 0 ) {
    $patch = 1;
    my $slices = $sa->fetch_by_region_unique( $slice->coord_system_name(),$slice->seq_region_name(),undef,undef,undef,$slice->coord_system()->version() );
    foreach my $slice ( @{$slices} ) {
      push @$genes,@{$ga->fetch_all_by_Slice($slice)};
      #      my $start = $slice->start;
    }
  }
  else {
    $genes = $ga->fetch_all_by_Slice($slice);
  }
  return ($genes, $patch);
}



=head2 get_attrib_values

  Arg[1]      : Arrayref of B::E::Attributes
  Arg[2]      : 'code' to search for
  Arg[3]      : 'value' to search for (optional)
  Example     : my $c = $self->get_attrib_values($attribs,'name'));
  Description : (i) In the absence of an attribute value argument, examines an arrayref
                of B::E::Attributes for a particular attribute type, returning the values
                for each attribute of that type. Can therefore be used to test for the
                number of attributes of that type.
                (ii) In the presence of the optional value argument it returns all
                attributes with that value ie can be used to test for the presence of an
                attribute with that particular value.
  Return type : arrayref of values for that attribute
  Caller      : general
  Status      : stable

=cut

sub get_attrib_values {
  my $self    = shift;
  my $attribs = shift;
  my $code    = shift;
  my $value   = shift;
  if (my @atts = grep {$_->code eq $code } @$attribs) {
    my $r = [];
    if ($value) {
      if (my @values = grep {$_->value eq $value} @atts) {
	foreach (@values) {
	  push @$r, $_->value;
	}
	return $r;
      }
      else {
	return [];
      }
    }
    else {
      foreach (@atts) {
	push @$r, $_->value;
      }
      return $r;
    }
  }
  else {
    return [];
  }
}

=head2 fix_attrib_value

  Arg[1]      : Arrayref of existing B::E::Attributes
  Arg[2]      : dbID of object
  Arg[3]      : name of object (just for reporting)
  Arg[4]      : attrib_type.code
  Arg[5]      : attrib_type.value
  Arg[6]      : interactive ? (0 by default)
  Arg[7]      : table
  Example     : $support->fix_attrib_value($attribs,$chr_id,$chr_name,'vega_export_mod','N',1);
  Description : adds a new attribute to an object, or updates an existing attribute with a new value
                Can be run in interactive or non-interactive mode (default)
  Return type : arrayref of results
  Caller      : general
  Status      : only ever tested with seq_region_attributes to date

=cut

sub fix_attrib_value {
  my $self        = shift;
  my $attribs     = shift;
  my $id          = shift;
  my $name        = shift;
  my $code        = shift;
  my $value       = shift;
  my $interact    = shift || 0;
  my $table       = shift || 'seq_region_attrib';

  #transiently set interactive parameter to zero
  my $int_before;
  if (! $interact) {
    $int_before = $self->param('interactive');
    $self->param('interactive',0);
  }

  #get any existing value(s) for this attribute
  my $existings = $self->get_attrib_values($attribs,$code);
	
  #add a new attribute if there is none...
  if (! @$existings ) {
    if ($self->user_proceed("Do you want to set $name attrib (code = $code) to value $value ?")) {
      my $r = $self->store_new_attribute($id,$code,$value);

      #reset interactive parameter
      $self->param('interactive',$int_before) if (! $interact);
      return $r;
		}
  }
  #...warn and exit if you're trying to update more than one value for the same attribute...
  elsif (scalar @$existings > 1) {
    $self->log_warning("You shouldn't be trying to update multiple attributes with the same code at once ($name:$code,$value), looks like you have duplicate entries in the (seq_region_)attrib table\n");
    exit;
  }

  #...or update an attribute with new values...
  else {
    my $existing = $existings->[0];
    if ($existing ne $value) {
      if ($self->user_proceed("Do you want to reset $name attrib (code = $code) from $existing to $value ?")) {
	my $r = $self->update_attribute($id,$code,$value);
	$self->param('interactive',$int_before) if (! $interact);
	push @$r, $existing;
	return $r;
      }
    }
    #...or make no change
    else {
      $self->param('interactive',$int_before) if (! $interact);
      return [];
    }
  }
}

=head2 _get_attrib_id

  Arg[1]      : attrib_type.code 
  Arg[2]      : database handle 
  Example     : $self->_get_attrib_id('name',$dbh)
  Description : get attrib_type.attrib_type_id from a attrib_type.code
  Return type : attrib_type.attrib_type_id 
  Caller      : internal
  Status      : stable

=cut

sub _get_attrib_id {
  my $self        = shift;
  my $attrib_code = shift;
  my $dbh         = shift;
  my ($attrib_id) = $dbh->selectrow_array(
    qq(select attrib_type_id
           from attrib_type
           where code = ?),
    {},
    ($attrib_code)
  );
  if (! $attrib_id) {
    $self->log_warning("There is no attrib_type_id for code $attrib_code, please patch the attrib_table\n");
    exit;
  }
  else {
    return $attrib_id;
  }
}

=head2 store_new_attribute

  Arg[1]      : seq_region.seq_region_id
  Arg[2]      : attrib_type.code
  Arg[3]      : attrib_type.value
  ARG[4]      : table to update (seq_region_attribute by default)
  Example     : $support->store_new_attribute(23,name,5);
  Description : uses MySQL to store an entry (code and value) in an attribute table 
                (seq_region_attrib by default)
  Return type : array_ref
  Caller      : general
  Status      : stable

=cut

sub store_new_attribute {
  my $self         = shift;
  my $sr_id        = shift;
  my $attrib_code  = shift;
  my $attrib_value = shift || '';
  my $table        = shift || 'seq_region_attrib';

  #get database handle
  my $dbh = $self->get_dbconnection('loutre');
  #get attrib_type_id for this particular attribute
  my $attrib_id = $self->_get_attrib_id($attrib_code,$dbh);
  #store
  my $r = $dbh->do(
    qq(insert into $table
           values (?,?,?)),
    {},
    ($sr_id,$attrib_id,$attrib_value)
  );
  return ['Stored',$r];
}

=head2 update_attribute

  Arg[1]      : seq_region.seq_region_id
  Arg[2]      : attrib_type.code
  Arg[3]      : attrib_type.value
  ARG[4]      : table to update (seq_region_attribute by default)
  Example     : $support->update_attribute(23,name,5);
  Description : uses MySQL to update an attribute table (seq_region_attrib by default)
  Return type : array_ref
  Caller      : general
  Status      : stable

=cut

sub update_attribute {
  my $self = shift;
  my $sr_id = shift;
  my $attrib_code  = shift;
  my $attrib_value = shift;
  my $table        = shift || 'seq_region_attrib';
  my $dbh = $self->get_dbconnection('loutre');
  my $attrib_id = $self->_get_attrib_id($attrib_code,$dbh);
  #update
  my $r = $dbh->do(
    qq(update $table
           set value = ?
           where seq_region_id = $sr_id
           and attrib_type_id = $attrib_id),
    {},
    ($attrib_value)
  );
  return ['Updated',$r];
}


=head2 remove_duplicate_attribs

  Arg[1]      : db handle
  Arg[2]      : table
  Example     : $support->remove_duplicate_attribs($dbh,'gene');
  Description : uses MySQL to remove duplicate entries from an attribute table
  Return type : none
  Caller      : general
  Status      : stable

=cut

sub remove_duplicate_attribs {
  my $self  = shift;
  my $dbh   = shift;
  my $table = shift;
  $dbh->do(qq(create table nondup_${table}_attrib like ${table}_attrib));
  $dbh->do(qq(insert into nondup_${table}_attrib (select ${table}_id, attrib_type_id, value from ${table}_attrib group by ${table}_id, attrib_type_id, value)));
  $dbh->do(qq(delete from ${table}_attrib));
  $dbh->do(qq(insert into ${table}_attrib (select ${table}_id, attrib_type_id, value from nondup_${table}_attrib)));
  $dbh->do(qq(drop table nondup_${table}_attrib));
}

=head2 sav_seq

  Arg[1]      : string (sequence to save)
  Example     : $support->save_seq('ACGT')
  Description : creates a temporary file containing the sequence you give it
  Return type : string (filename)
  Caller      : general
  Status      : stable

=cut

sub save_seq {
  my $self = shift;
  my $content = shift ;
  my $seq_file = $self->param('logpath') . '/SEQ_' . time() . int(rand()*100000000) . $$;
  open (my $fh,">$seq_file", $seq_file) or die("Cannot create working file.$!");
  print $fh $content;
  close $fh;
  return ($seq_file);
}

=head2 get_alignment

  Arg[1]      : string (first sequence)
  Arg[1]      : string (second sequence)
  Arg[1]      : string (sequence type))
  Example     : $support->get_alignment('AAAAA','CCCCCCC','DNA')
  Description : creates a temporary file containing the sequence you give it
  Return type : string (filename)
  Caller      : general
  Status      : stable

=cut

sub get_alignment {
  my $self = shift;
  my $ext_seq  = shift || return undef;
  my $int_seq  = shift || return undef;
  $int_seq =~ s/<br \/>//g;
  my $seq_type = shift || return undef;


  # To stop box running out of memory - put an upper limit on the size of sequence
  # that alignview can handle
  if (length $int_seq > 1e6 || length $ext_seq > 1e6)  {
    $self->log_error('Cannot align if sequence > 1 Mbase');
  }

  my $int_seq_file = $self->save_seq($int_seq);
  my $ext_seq_file = $self->save_seq($ext_seq);

  my $label_width  = '22'; # width of column for e! object label
  my $output_width = 61;   # width of alignment
  my $dnaAlignExe  = '/localsw/bin/emboss/bin/matcher -asequence %s -bsequence %s -outfile %s';
  my $pepAlignExe  = '/localsw/bin/wise2/bin/psw -dymem explicit -m /localsw/bin/wise2/wisecfg/blosum62.bla %s %s -n %s -w %s > %s';

  my $out_file = time() . int(rand()*100000000) . $$;
  $out_file = $self->param('logpath').'/' . $out_file . '.out';

  my $command;
  my $fh;
  if ($seq_type eq 'DNA') {
    $command = sprintf $dnaAlignExe, $int_seq_file, $ext_seq_file, $out_file;
    `$command`;
    unless (open($fh, "<", $out_file)) {
      $command = sprintf $dnaAlignExe, $int_seq_file, $ext_seq_file, $out_file;
      `$command`;
    }
  }
  elsif ($seq_type eq 'PEP') {
    $command = sprintf $pepAlignExe, $int_seq_file, $ext_seq_file, $label_width, $output_width, $out_file;
    `$command`;
    unless (open($fh, "<", $out_file)) {
      $self->log_warning("Cannot open alignment file\n");
    }
  }
  my $alignment ;
  while (<$fh>) {
    next if $_ =~ 
    /\#Report_file
     |\#----.*
     |\/\/\s*
     |\#\#\#
     |^\#$
     |Rundate: #matcher
     |Commandline #matcher
     |asequence #matcher
     |bsequence #matcher
     |outfile #matcher
     |aformat #matcher
     |Align_format #matcher
     |Report_file #matcher
     /x;
    $alignment .= $_;
  }

  $alignment =~ s/\n+$//;
  unlink $out_file;
  unlink $int_seq_file;
  unlink $ext_seq_file;
  return $alignment;
}

sub allowed_duplicate_regions {
  my $self = shift;

  # set up lists of vega names of wanted haplotype / strains
  my $regions = {
    'human' => [
      {
        '6'      => '28000000:34000000',
        '6-COX'  => 'all',
        '6-QBL'  => 'all',
        '6-APD'  => 'all',
        '6-MANN' => 'all',
        '6-MCF'  => 'all',
        '6-DBB'  => 'all',
        '6-SSTO' => 'all',
      },
      {
        '19'       => '54020000:54910000',
        '19-PGF_1' => 'all',
        '19-PGF_2' => 'all',
        '19-COX_1' => 'all',
        '19-COX_2' => 'all',
        '19-DM1A'  => 'all',
        '19-DM1B'  => 'all',
        '19-MC1A'  => 'all',
        '19-MC1B'  => 'all',
      },
    ],
    'mouse' => [
      {
        '1'               => '60564732:63711641',
        '1-Idd5.1_DIL'    => 'all',
        '1-Idd5.1_CHO'    => 'all',
      },
      {
        '1'               => '65533102:69307244',
        '1-Idd5.3_DIL'    => 'all',
      },
      {
        '1'               => '130232728:130661594',
        '1-Idd5.4_DIL'    => 'all',
      },
      {
        '3'               => '36492618:37600833',
        '3-Idd3_DIL'      => 'all',
        '3-Idd3_129'      => 'all',
      },
      {
        '3'               => '99848826:101467080',
        '3-Idd10_DIL'     => 'all',
      },
      {
        '3'               => '109144756:109930492',
        '3-Idd18.1_DIL'   => 'all',
      },
      {
        '3'               => '103489414:104054885',
        '3-Idd18.2_DIL'   => 'all',
      },
      {
        '4'               => '128371876:131853368',
        '4-Idd9.1_DIL'    => 'all',
      },
      {
        '4'               => '134841437:135252443',
        '4-Idd9.1M_DIL'   => 'all',
      },
      {
        '4'               => '146049976:149895141',
        '4-Idd9.2_DIL'    => 'all',
      },
      {
        '4'               => '149556939:151385803',
        '4-Idd9.3_DIL'    => 'all',
      },
      {
        '6'               => '143550839:149565172',
        '6-Idd6.1_2_CHO'  => 'all',
      },
      {
        '6'               => '129593784:131241919',
        '6-Idd6.AM_CHO'   => 'all',
      },
      {
        '11'              => '69890356:71340035',
        '11-Idd4.1_DIL'   => 'all'
      },
      {
        '11'              => '72734492:74404570',
        '11-Idd4.2_DIL'   => 'all',
      },	
      {
        '11'              => '86785996:90007691',
        '11-Idd4.2Q_CHO'  => 'all',
      },
      {
        '17'              => '27302611:29220265',
        '17-Idd16.1_CHO'  => 'all',
      },
      {
        '17'              => '33567721:38721962',
        '17-Idd1_CHO'     => 'all',
        '17-Idd1_DIL'     => 'all',
      },
    ],
    'pig' => [
      {
        '7'               => '24728583:29807435',
        '7-LW'            => 'all',
      },
      {
        'X'               => 'all',
        'X-WTSI'          => 'all',
      },
    ],
    'zebrafish' => [],
    'gorilla' => [],
    'wallaby' => [],
    'chimp' => [],
    'rat'   => [],
  };

return $regions;

}

1;

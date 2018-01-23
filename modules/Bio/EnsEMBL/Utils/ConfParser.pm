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

Bio::EnsEMBL::Utils::ConfParser - configuration parser for perl scripts

=head1 SYNOPSIS

  my $conf = new Bio::EnsEMBL::Utils::ConfParser(
    -SERVERROOT   => "/path/to/ensembl",
    -DEFAULT_CONF => "my.default.conf"
  );

  # parse options from configuration file and commandline
  $conf->parse_options(
    'mandatory_string_opt=s' => 1,
    'optional_numeric_opt=n' => 0,
  );

  # get a paramter value
  my $val = $conf->param('manadatory_string_op');

=head1 DESCRIPTION

This module parses a configuration file and the commandline options
passed to a script (the latter superseed the former). Configuration
files contain ini-file style name-value pairs, and the commandline
options are passed to Getopt::Long for parsing.

The parameter values are consequently accessible via the param()
method. You can also create a commandline string of all current
parameters and their values to pass to another script.

=cut

package Bio::EnsEMBL::Utils::ConfParser;

use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use Text::Wrap;
use Cwd qw(abs_path);
use Pod::Usage qw(pod2usage);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(user_proceed);


=head2 new

  Arg [SERVERROOT] :
                String $serverroot - root directory of your ensembl code
  Arg [DEFAULT_CONF] :
                String $default_conf - default configuration file
  Example     : my $conf = new Bio::EnsEMBL::Utils::ConfParser(
                  -SERVERROOT => '/path/to/ensembl',
                  -DEFAULT_CONF => 'my.default.conf'
                );
  Description : object constructor
  Return type : Bio::EnsEMBL::Utils::ConfParser object
  Exceptions  : thrown if no serverroot is provided
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($serverroot, $default_conf) =
    rearrange([qw(SERVERROOT DEFAULT_CONF)], @_);

  throw("You must supply a serverroot.") unless ($serverroot);

  my $self = {};
  bless ($self, $class);

  $self->serverroot($serverroot);
  $self->default_conf($default_conf || "$ENV{HOME}/.ensembl_script.conf");

  return $self;
}


=head2 parse_options

  Arg[1..n]   : pairs of option definitions and mandatory flag (see below for
                details)
  Example     : $conf->parse_options(
                  'mandatory_string_opt=s' => 1,
                  'optional_numeric_opt=n' => 0,
                );
  Description : This method reads options from an (optional) configuration file
                and parses the commandline options supplied by the user.
                Commandline options will superseed config file settings. The
                string "$SERVERROOT" in the configuration entries will be
                replaced by  the appropriate value.

                The arguments passed to this method are pairs of a Getopt::Long
                style option definition (in fact it will be passed to
                GetOptions() directly) and a flag indicating whether this
                option is mandatory (1) or optional (0).

                In addition to these user-defined options, a set of common
                options is always parsed. See _common_options() for details.
                
                If you run your script with --interactive the user will be
                asked to confirm the parameters after parsing.
                
                All parameters will then be accessible via $self->param('name').
  Return type : true on success 
  Exceptions  : thrown if configuration file can't be opened
                thrown on missing mandatory parameters
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub parse_options {
  my ($self, @params) = @_;

  # add common options to user supplied list
  push @params, $self->_common_options;

  # read common commandline options
  my %h;
  my %params = @params;

  Getopt::Long::Configure('pass_through');
  &GetOptions(\%h, keys %params);

  # reads config file
  my $conffile = $h{'conffile'} || $self->default_conf;
  $conffile = abs_path($conffile);

  if (-e $conffile) {
    open(my $fh, '<', $conffile) or throw( 
        "Unable to open configuration file $conffile for reading: $!");

    my $serverroot = $self->serverroot;
    my $last;

    while (my $line = <$fh>) {
      chomp $line;
      
      # remove leading and trailing whitespace
      $line =~ s/^\s*//;
      $line =~ s/\s*$//;

      # join with next line if terminated with backslash (this is to allow
      # multiline configuration settings
      $line = $last . $line;
      if ($line =~ /\\$/) {
        $line =~ s/\\$//;
        $last = $line;
        next;
      } else {
        $last = undef;
      }

      # remove comments
      $line =~ s/^[#;].*//;
      $line =~ s/\s+[;].*$//;

      # read options into internal parameter datastructure
      next unless ($line =~ /(\w\S*)\s*=\s*(.*)/);
      my $name = $1;
      my $val = $2;

      # strip optional quotes from parameter values
      $val =~ s/^["'](.*)["']/$1/;

      # replace $SERVERROOT with value
      if ($val =~ /\$SERVERROOT/) {
        $val =~ s/\$SERVERROOT/$serverroot/g;
        $val = abs_path($val);
      }
      $self->param($name, $val);
    }
    close($fh);

    $self->param('conffile', $conffile);
  }

  # override configured parameter with commandline options
  map { $self->param($_, $h{$_}) } keys %h;

  # check for required params, convert comma to list, maintain an ordered
  # list of parameters and list of 'flag' type params
  my @missing = ();
  my $i = 0;

  foreach my $param (@params) {
    next if ($i++ % 2);

    my $required = $params{$param};
    my ($list, $flag);
    $list = 1 if ($param =~ /\@$/);
    $flag = 1 if ($param =~ /!$/);
    $param =~ s/(^\w+).*/$1/;
    
    $self->comma_to_list($param) if ($list);

    push @missing, $param if ($required and !$self->param($param));
    push @{ $self->{'_ordered_params'} }, $param;
    $self->{'_flag_params'}->{$param} = 1 if ($flag);
  }
  
  if (@missing) {
    throw("Missing parameters: @missing.\nYou must specify them on the commandline or in your conffile.\n");
  }

  # error handling and --help
  pod2usage(1) if ($self->param('help'));

  # ask user to confirm parameters to proceed
  $self->confirm_params;

  return(1);
}


#
# Commonly used options. These are parsed by default even if they are not
# passed to parse_options() explicitely.
#
sub _common_options {
  my $self = shift;
  return (
    'conffile|conf=s' => 0,
    'logfile|log=s' => 0,
    'logauto!' => 0,
    'logautobase=s' => 0,
    'logautoid=s' => 0,
    'logpath=s' => 0,
    'logappend|log_append|log-append!' => 0,
    'loglevel=s' => 0,
    'is_component|is-component!' => 0,
    'interactive|i!' => 0,
    'dry_run|dry-run|dry|n!' => 0,
    'help|h|?' => 0,
  );
}


=head2 confirm_params

  Example     : $conf->confirm_params;
  Description : If the script is run with the --interactive switch, this method
                prints a table of all parameters and their values and asks user
                to confirm if he wants to proceed.
  Return type : true on success
  Exceptions  : none
  Caller      : parse_options()
  Status      : At Risk
              : under development

=cut

sub confirm_params {
  my $self = shift;

  if ($self->param('interactive')) {
    # print parameter table
    print "Running script with these parameters:\n\n";
    print $self->list_param_values;

    # ask user if he wants to proceed
    exit unless user_proceed("Continue?", 1, 'n');
  }
  
  return(1);
}


=head2 param

  Arg[1]      : Parameter name
  Arg[2..n]   : (optional) List of values to set
  Example     : # getter
                my $dbname = $conf->param('dbname');

                # setter
                $conf->param('port', 3306);
                $conf->param('chromosomes', 1, 6, 'X');
  Description : Getter/setter for parameters. Accepts single-value params and
                list params.
  Return type : Scalar value for single-value parameters, array of values for
                list parameters
  Exceptions  : thrown if no parameter name is supplied
  Caller      : general
  Status      : At Risk
              : under development

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


=head2 is_true

  Arg[1]      : Parameter name
  Example     : unless ($conf->is_true('upload')) {
                  print "Won't upload data.\n";
                  next;
                }
  Description : Checks whether a param value is set to 'true', which is defined
                here as TRUE (in the Perl sense) but not the string 'no'.
  Return type : Boolean
  Exceptions  : thrown if no parameter name is supplied
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_true {
  my $self = shift;
  my $name = shift or throw("You must supply a parameter name");

  my $param = $self->param($name);

  if ($param and !($param =~ /^no$/i)) {
    return(1);
  } else {
    return(0);
  }
}


=head2 list_params

  Example     : print "Current parameter names:\n";
                foreach my $param (@{ $conf->list_params }) {
                  print "  $param\n";
                }
  Description : Returns a list of the currently available parameter names. The
                list will be in the same order as option definitions were
                passed to the new() method.
  Return type : Arrayref of parameter names
  Exceptions  : none
  Caller      : list_param_values(), create_commandline_options()
  Status      : At Risk
              : under development

=cut

sub list_params {
  my $self = shift;
  return $self->{'_ordered_params'} || [];
}


=head2 list_param_values

  Example     : print LOG $conf->list_param_values;
  Description : prints a table of the parameters used in the script
  Return type : String - the table to print
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub list_param_values {
  my $self = shift;

  $Text::Wrap::colums = 72;

  my $txt = sprintf "    %-20s%-40s\n", qw(PARAMETER VALUE);
  $txt .= "    " . "-"x70 . "\n";

  foreach my $key (@{ $self->list_params }) {
    my $val;
    if (defined($self->param($key))) {
      $txt .= Text::Wrap::wrap(sprintf('    %-19s ', $key), ' 'x24,
        join(", ", $self->param($key)))."\n";
    }
  }

  $txt .= "\n";

  return $txt;
}


=head2 create_commandline_options

  Arg[1..n]   : param/value pairs which should be added to or override the
                currently defined parameters
  Example     : $conf->create_commandline_options(
                    'dbname' => 'homo_sapiens_vega_33_35e',
                    'interactive' => 0
                );
  Description : Creates a commandline options string of all current paramters
                that can be passed to another script.
  Return type : String - commandline options string
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub create_commandline_options {
  my ($self, %replace) = @_;

  my %param_hash;

  # deal with list values
  foreach my $param (@{ $self->list_params }) {
    my ($first, @rest) = $self->param($param);
    next unless (defined($first));

    if (@rest) {
      $first = join(",", $first, @rest);
    }
    $param_hash{$param} = $first;
  }

  # replace values
  foreach my $key (keys %replace) {
    $param_hash{$key} = $replace{$key};
  }

  # create the commandline options string
  my $options_string;
  foreach my $param (keys %param_hash) {

    my $val = $param_hash{$param};

    # deal with 'flag' type params correctly
    if ($self->{'_flag_params'}->{$param}) {
      # change 'myparam' to 'nomyparam' if no value set
      $param = 'no'.$param unless ($val);

      # unset value (this is how flags behave)
      $val = undef;
    } else {
      # don't add the param if it's not a flag param and no value is set
      next unless (defined($val));

      # quote the value if it contains blanks
      if ($val =~ /\s+/) {
        # use an appropriate quoting style
        ($val =~ /'/) ? ($val = qq("$val")) : ($val = qq('$val'));
      }
    }
    
    $options_string .= sprintf(qq(--%s %s ), $param, $val);
  }

  return $options_string;
}


=head2 comma_to_list

  Arg[1..n]   : list of parameter names to parse
  Example     : $conf->comma_to_list('chromosomes');
  Description : Transparently converts comma-separated lists into arrays (to
                allow different styles of commandline options, see perldoc
                Getopt::Long for details). Parameters are converted in place
                (accessible through $self->param('name')).
  Return type : true on success
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

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
  Example     : $conf->list_or_file('gene');
  Description : Determines whether a parameter holds a list or it is a filename
                to read the list entries from.
  Return type : true on success
  Exceptions  : thrown if list file can't be opened
  Caller      : general
  Status      : At Risk
              : under development

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


=head2 serverroot

  Arg[1]      : (optional) String - root directory of your ensembl checkout
  Example     : my $serverroot = $conf->serverroot;
  Description : Getter/setter for the root directory of your ensembl checkout.
  Return type : String
  Exceptions  : none
  Caller      : new(), general
  Status      : At Risk
              : under development

=cut

sub serverroot {
  my $self = shift;
  $self->{'_serverroot'} = shift if (@_);
  return $self->{'_serverroot'};
}


=head2 default_conf

  Arg[1]      : (optional) String - default configuration file
  Example     : $conf->default_conf('my.default.conf');
  Description : Getter/setter for the default configuration file.
  Return type : String
  Exceptions  : none
  Caller      : new(), general
  Status      : At Risk
              : under development

=cut

sub default_conf {
  my $self = shift;
  $self->{'_default_conf'} = shift if (@_);
  return $self->{'_default_conf'};
}


1;


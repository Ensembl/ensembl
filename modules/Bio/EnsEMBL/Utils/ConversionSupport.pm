package Bio::EnsEMBL::Utils::ConversionSupport;

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
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use FindBin qw($Bin $Script);
use POSIX qw(strftime);

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
    (my $serverroot = shift) or throw("You must supply a serverroot");
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
        'driver|dbdriver|db_driver=s',
        'conffile|conf=s',
        'logfile|log=s',
        'logpath=s',
        'interactive|i=s',
        'dry_run|dry|n=s',
        'help|h|?',
    );

    # reads config file
    my $conffile = $h{'conffile'} || $self->serverroot . "/conf/Conversion.ini";
    if (-e $conffile) {
        open(CONF, $conffile) or throw( 
            "Unable to open configuration file $conffile for reading: $!");
        while (<CONF>) {
            chomp;

            # remove comments
            s/^[#;].*//;
            s/\s+[;].*$//;

            # read options into internal parameter datastructure
            next unless (/(\w\S*)\s*=\s*(.*)/);
            $self->param($1, $2);
        }
    } else {
        warning("Unable to open configuration file $conffile for reading: $!");
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
    $self->error($@) if $@;
    return(1);
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
    print "Running script with these parameters:\n\n";
    print $self->list_all_params;

    # ask user if he wants to proceed
    $self->user_confirm;
    
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
    my $txt = sprintf "    %-20s%-40s\n", qw(PARAMETER VALUE);
    $txt .= "    " . "-"x70 . "\n";
    $Text::Wrap::colums = 72;
    foreach my $key (sort keys %{ $self->{'_param'} }) {
        my @vals = $self->param($key);
        $txt .= Text::Wrap::wrap( sprintf('    %-20s', $key),
                                  ' 'x24,
                                  join(", ", @vals)
                                ) . "\n";
    }
    $txt .= "\n";
    return $txt;
}

=head2 user_confirm

  Example     : print "Do you want to continue?\n";
                $support->user_confirm;
  Description : If running interactively, the user is asked if he wants to
                proceed.
  Return type : true on success.
  Exceptions  : none
  Caller      : general

=cut

sub user_confirm {
    my $self = shift;

    if ($self->param('interactive')) {
        print "Continue? [y/N] ";
        my $input = lc(<>);
        chomp $input;
        unless ($input eq 'y') {
            print "Aborting.\n";
            exit(0);
        }
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
        $self->param($param,
            split (/,/, join (',', $self->param($param))));
    }
    return(1);
}

=head2 list_to_file

  Arg[1]      : Name of parameter to parse
  Example     : $support->list_to_file('gene_stable_id');
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
  Example     : my $db = $support->get_database('core');
  Description : Connects to the database specified.
  Return type : DBAdaptor of the appropriate type
  Exceptions  : thrown if asking for unknown database
  Caller      : general

=cut

sub get_database {
    my $self = shift;
    my $database = shift or throw("You must provide a database");
    $self->check_required_params(qw(host port user pass dbname));

    my %adaptors = (
        core    => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
        ensembl => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
        otter   => 'Bio::Otter::DBSQL::DBAdaptor',
        vega    => 'Bio::Otter::DBSQL::DBAdaptor',
    );
    my %valid = map { $_ => 1 } keys %adaptors;
    throw("Unknown database: $database") unless $valid{$database};

    $self->dynamic_use($adaptors{$database});
    my $dba = $adaptors{$database}->new(
            -host   => $self->param('host'),
            -port   => $self->param('port'),
            -user   => $self->param('user'),
            -pass   => $self->param('pass'),
            -dbname => $self->param('dbname'),
    );
    return $dba;
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
    my ($parent_namespace, $module) = $classname =~/^(.*::)(.*)$/ ?
                                        ($1,$2) : ('::', $classname);
    no strict 'refs';
    # return if module has already been imported
    return 1 if $parent_namespace->{$module.'::'};
    
    eval "require $classname";
    throw("Failed to require $classname: $@") if ($@);
    $classname->import();
    
    return 1;
}

=head2 get_chrlength

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBAdaptor $dba
  Example     : my $chr_length = $support->get_chrlength($dba);
  Description : Get all chromosomes and their length from the database. Return
                chr_name/length for the chromosomes the user requested (or all
                chromosomes by default)
  Return type : Hashref - chromosome_name => length
  Exceptions  : thrown if not passing a Bio::EnsEMBL::DBSQL::DBAdaptor
  Caller      : general

=cut

sub get_chrlength {
    my ($self, $dba) = @_;
    throw("get_chrlength should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n")
        unless ($dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));

    my $sa = $dba->get_SliceAdaptor;
    my @chromosomes = map { $_->seq_region_name } 
                            @{ $sa->fetch_all('chromosome') };
    my %chr = map { $_ => $sa->fetch_by_region('chromosome', $_)->length }
                            @chromosomes;

    my @wanted = $self->param('chromosomes');
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
  Exceptions  : thrown if species name can't be determined from db
  Caller      : general

=cut

sub get_species_scientific_name {
    my ($self, $dba) = @_;
    my $sql = qq(
        SELECT
                meta_value
        FROM
                meta
        WHERE meta_key = "species.classification"
        ORDER BY meta_id
        LIMIT 2
    );
    my $sth = $dba->dbc->db_handle->prepare($sql);
    $sth->execute;
    my @sp;
    while (my @row = $sth->fetchrow_array) {
        push @sp, $row[0];
    }
    $sth->finish;
    my $species = join(" ", reverse @sp);
    $self->throw("Could not determine species scientific name from database.")
        unless $species;
    return $species;
}

=head2 sort_chromosomes

  Arg[1]      : Hashref $chr_hashref - Hashref with chr_name as keys
  Example     : my $chr = { '6-COX' => 1, '1' => 1, 'X' => 1 };
                my @sorted = $support->sort_chromosomes($chr);
  Description : Sorts chromosomes in an intuitive way (numerically, then
                alphabetically)
  Return type : List - sorted chromosome names
  Exceptions  : thrown if no hashref is provided
  Caller      : general

=cut

sub sort_chromosomes {
    my ($self, $chr_hashref) = @_;
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

=head2 log

  Arg[1]      : String $txt - the text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle('>>');
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
    $txt = "    "x$indent . $txt;
    my $fh = $self->{'_log_filehandle'};
    throw("Unable to obtain log filehandle") unless $fh;
    print $fh "$txt";
    return(1);
}

=head2 log_warning

  Arg[1]      : String $txt - the warning text to log
  Arg[2]      : Int $indent - indentation level for log message
  Example     : my $log = $support->log_filehandle('>>');
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

=head2 log_filehandle

  Arg[1]      : String $mode - file access mode
  Example     : my $log = $support->log_filehandle('>>');
                # print to the filehandle
                print $log 'Lets start logging...\n';
                # log via the wrapper $self->log()
                $support->log('Another log message.\n');
  Description : Returns a filehandle for logging (STDERR by default, logfile if
                set from config or commandline). You can use the filehandle
                directly to print to, or use the smart wrapper $self->log()
  Return type : Filehandle - the filehandle to log to
  Exceptions  : thrown if logfile can't be opened
  Caller      : general

=cut

sub log_filehandle {
    my ($self, $mode) = @_;
    $mode ||= ">";
    my $fh = \*STDERR;
    if (my $logfile = $self->param('logfile')) {
        if (my $logpath = $self->param('logpath')) {
            $logfile = "$logpath/$logfile";
        }
        open($fh, "$mode", $logfile) or throw(
            "Unable to open $logfile for writing: $!");
    }
    $self->{'_log_filehandle'} = $fh;
    return $self->{'_log_filehandle'};
}

=head2 init_log

  Example     : print LOG $support->init_log;
  Description : Returns some header information for a logfile. This includes
                script name, date, user running the script and parameters the
                script will be running with
  Return type : String - the log text
  Exceptions  : none
  Caller      : general

=cut

sub init_log {
    my $self = shift;

    # print script name, date, user who is running it
    my $hostname = `hostname`;
    chomp $hostname;
    my $script = "$hostname:$Bin/$Script";
    my $user = `whoami`;
    chomp $user;
    my $txt = "Script: $script\nDate: ".$self->date."\nUser: $user\n";

    # print parameters the script is running with
    $txt .= "Parameters:\n\n";
    $txt .= $self->list_all_params;

    return $txt;
}

=head2 finish_log

  Example     : print LOG $support->finish_log;
  Description : Return footer information to write to a logfile. This includes
                the number of logged warnings, timestamp and memory footprint.
  Return type : String - the log text
  Exceptions  : none
  Caller      : general

=cut

sub finish_log {
    my $self = shift;
    my $txt = "All done. ".$self->warnings." warnings. ".$self->date_and_mem."\n";
    return $txt;
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

1;

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

Bio::EnsEMBL::Utils::CliHelper

=head1 VERSION

$Revision$

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::CliHelper;

  my $cli = Bio::EnsEMBL::Utils::CliHelper->new();

  # get the basic options for connecting to a database server
  my $optsd = $cli->get_dba_opts();

  # add another option
  push(@$optsd,"print");

  # process the command line with the supplied options plus a reference to a help subroutine
  my $opts = $cli->process_args($optsd,\&usage);
  
  # use the command line options to get an array of database details
  for my $db_args (@{$cli->get_dba_args_for_opts($opts)}) {
    # use the args to create a DBA
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{$db_args});
    ...
    if(defined $opts->{print}) {
    	...
    }
  }
  
  For adding secondary databases, a prefix can be supplied. For instance, to add a second set of
  db params prefixed with dna (-dnahost -dbport etc.) use the prefix argument with get_dba_opts and 
  get_dba_args_for_opts:
  # get the basic options for connecting to a database server
  my $optsd =
   [ @{ $cli_helper->get_dba_opts() }, @{ $cli_helper->get_dba_opts('gc') } ];
  # process the command line with the supplied options plus a help subroutine
  my $opts = $cli_helper->process_args( $optsd, \&usage );
  # get the dna details
  my ($dna_dba_details) =
    @{ $cli_helper->get_dba_args_for_opts( $opts, 1, 'dna' ) };
  my $dna_db =
    Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{$dna_dba_details} ) );

=head1 DESCRIPTION

Utilities for a more consistent approach to parsing and handling EnsEMBL script command lines

=head1 METHODS

See subroutines.

=cut

package Bio::EnsEMBL::Utils::CliHelper;

use warnings;
use strict;

use Carp;
use Getopt::Long qw(:config auto_version no_ignore_case);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $dba_opts = [{args => ['host', 'dbhost', 'h'], type => '=s'},
				{args => ['port', 'dbport', 'P'], type => ':i'},
				{args => ['user', 'dbuser', 'u'], type => '=s'},
				{args => ['pass', 'dbpass', 'p'], type => ':s'},
				{args => ['dbname',  'D'],         type => ':s'},
				{args => ['pattern', 'dbpattern'], type => ':s'},
				{args => ['driver'],     type => ':s'},
				{args => ['species_id'], type => ':i'},
				{args => ['species'],    type => ':i'},];

=head2 new()

  Description : Construct a new instance of a CliHelper object
  Returntype  : Bio::EnsEMBL::Utils:CliHelper
  Status      : Under development

=cut

sub new {
  my ($class, @args) = @_;
  my $self = bless({}, ref($class) || $class);
  return $self;
}

=head2 get_dba_opts()

  Arg [1]     : Optional prefix for dbnames e.g. dna
  Description : Retrieves the standard options for connecting to one or more Ensembl databases
  Returntype  : Arrayref of option definitions
  Status      : Under development

=cut

sub get_dba_opts {
  my ($self, $prefix) = @_;
  $prefix ||= '';
  my @dba_opts = map {
	my $opt = join '|', map { $prefix . $_ } @{$_->{args}};
	$opt . $_->{type};
  } @{$dba_opts};
  return \@dba_opts;
}

=head2 process_args()

    Arg [1]     : Arrayref of supported command line options (e.g. from get_dba_opts)
    Arg [2]     : Ref to subroutine to be invoked when -help or -? is supplied
    Description : Retrieves the standard options for connecting to one or more Ensembl databases
    Returntype  : Hashref of parsed options
    Status      : Under development

=cut

sub process_args {
  my ($self, $opts_def, $usage_sub) = @_;
  my $opts = {};
  push @{$opts_def}, q/help|?/ => $usage_sub;
  GetOptions($opts, @{$opts_def}) ||
	croak 'Could not parse command line arguments';
  return $opts;
}

=head2 get_dba_args_for_opts()

    Arg [1]     : Hash of options (e.g. parsed from command line options by process_args())
    Arg [2]     : If set to 1, the databases are assumed to have a single species only. Default is 0 if database name matches collection, 1 otherwise.
    Arg [3]     : Optional prefix to use when parsing e.g. dna
    Description : Uses the parsed command line options to generate an array of DBAdaptor arguments 
                : (e.g. expands dbpattern, finds all species_ids for multispecies databases)
                : These can then be passed directly to Bio::EnsEMBL::DBSQL::DBAdaptor->new()
    Returntype  : Arrayref of DBA argument hash refs 
    Status      : Under development

=cut

sub get_dba_args_for_opts {
  my ($self, $opts, $single_species_opt, $prefix) = @_;
  $prefix ||= '';

  my ($host,    $port,   $user,    $pass, $dbname,
	  $pattern, $driver, $species, $species_id)
	= map { $prefix . $_ }
	qw(host port user pass dbname pattern driver species species_id);
	
  my @db_args;
  if (defined $opts->{$host}) {
	my $dbc =
	  Bio::EnsEMBL::DBSQL::DBConnection->new(-USER   => $opts->{$user},
											 -PASS   => $opts->{$pass},
											 -HOST   => $opts->{$host},
											 -PORT   => $opts->{$port},
											 -DRIVER => $opts->{$driver}
	  );
	my @dbnames;
	if (defined $opts->{$dbname}) {
	  push @dbnames, $opts->{$dbname};
	}
	elsif (defined $opts->{$pattern}) {
   # get a basic DBConnection and use to find out which dbs are involved
	  @dbnames = grep { m/$opts->{$pattern}/smx }
		@{$dbc->sql_helper()->execute_simple(q/SHOW DATABASES/)};
	}
	else {
	  croak 'dbname or dbpattern arguments required';
	}
	for my $dbname (@dbnames) {

#Decipher group of DBAdaptor by capturing the name_name(_name?)_core_ code. Otherwise we don't know
	  my ($group) = $dbname =~
		/^[a-z]+_[a-z0-9]+(?:_[a-z0-9]+)?_([a-z]+)(?:_\d+)?_\d+/;
	  # set multi where we have collections
	  my $multi = $dbname =~ m/_collection_/ ? 1 : 0;
	  my $species_ids;
	  my $single_species = $single_species_opt;

	  if (!defined $single_species) {
        # if we're dealing with a collection, turn off single species mode by default
		$single_species = $dbname =~ m/_collection_/ ? 0 : 1;
	  }
	  if ($single_species != 1) {
	  	# for multispecies, get the list of species from meta
		$species_ids =
		  $dbc->sql_helper()
		  ->execute(
"SELECT species_id,meta_value FROM $dbname.meta WHERE meta_key='species.production_name'"
		  );
		if (!defined $opts->{$species_id} &&
			scalar(@{$species_ids}) == 0)
		{
		  croak "No species.production_name found in database";
		}
	  }
	  # if we didn't get a list from meta, go ahead and use the supplied arguments if we have them
	  if (defined $opts->{$species_id}) {
		$species_ids = [[$opts->{$species_id}, $opts->{$species}]];
	  }
	  # otherwise assume the default species
	  elsif(!defined $species_ids) {
		$species_ids = [[1, undef]];
	  }
	  # deal with each species in turn
	  for my $species_id (@{$species_ids}) {
		my $args = {-HOST            => $opts->{$host},
					-USER            => $opts->{$user},
					-PORT            => $opts->{$port},
					-PASS            => $opts->{$pass},
					-DBNAME          => $dbname,
					-DRIVER          => $opts->{$driver},
					-SPECIES_ID      => $species_id->[0],
					-SPECIES         => $species_id->[1],
					-MULTISPECIES_DB => $multi};
		$args->{-GROUP} = $group if $group;
		push(@db_args, $args);
	  }
	} ## end for my $dbname (@dbnames)
  } ## end if (defined $opts->{$host...})
  else {
	croak '(db)host arguments required';
  }
  return \@db_args;
} ## end sub get_dba_args_for_opts

=head2 get_dba_args_for_opts()

    Arg [1]     : Hash of options (e.g. parsed from command line options by process_args())
    Arg [2]     : If set to 1, the databases are assumed to have a single species only. Default is 0.
    Arg [3]     : Optional prefix to use when parsing e.g. dna
    Description : Uses the parsed command line options to generate an array DBAdaptors. 
                : Note this can overload connections on a server
    Returntype  : Arrayref of Bio::EnsEMBL::DBSQL::DBAdaptor
    Status      : Under development

=cut

sub get_dbas_for_opts {
  my ($self, $opts, $single_species, $prefix) = @_;

# get all the DBA details that we want to work with and create DBAs for each in turn
  my $dbas;
  for my $args (
	   @{$self->get_dba_args_for_opts($opts, $single_species, $prefix)})
  {
	push @{$dbas}, Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$args});
  }
  return $dbas;
}

=head2 load_registry_for_opts

  Arg [1]    	: Hash of options (e.g. parsed from command line options by process_args()) 
  Arg [2]     : Optional prefix to use when parsing e.g. dna or master 
  Description	: Loads a Registry from the given options hash. If a C<registry> 
                option is given then the code will call C<load_all>. Otherwise
                we use the database parameters given to call 
                C<load_registry_from_db()>.
  Returntype  : Integer of the number of DBAdaptors loaded
  Status      : Under development

=cut

sub load_registry_for_opts {
  my ($self, $opts, $prefix) = @_;
  $prefix ||= q{};
  if ($opts->{registry}) {
	my $location = $opts->{registry};
	return Bio::EnsEMBL::Registry->load_all($location);
  }
  my ($host, $port, $user, $pass) =
	map { $prefix . $_ } qw(host port user pass);
  my %args = (-HOST => $opts->{$host},
			  -PORT => $opts->{$port},
			  -USER => $opts->{$user},);
  $args{-PASS} = $opts->{$pass};
  return Bio::EnsEMBL::Registry->load_registry_from_db(%args);
}

1;

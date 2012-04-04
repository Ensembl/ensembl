#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw( :config no_ignore_case );
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->defaults();
  my $databases = $self->databases();
  foreach my $db (@{$databases}) {
    $self->process_db($db);
  }
  return;
}

sub args {
  my ($self) = @_;
  my $opts = {};
  my @cmd_opts = qw/
    mhost|mh=s
    mport|mP=i
    muser|mu=s
    mpass|mp=s
    mdatabase|md=s
    host|h=s
    port|P=i
    user|u=s
    pass|p=s
    database|d=s
    pattern=s
    retirealiases!
    verbose|v!
    help
    man
  /;
  GetOptions($opts, @cmd_opts) or pod2usage(-verbose => 1, -exitval => 1);
  pod2usage(-verbose => 1, -exitval => 0) if $opts->{help};
  pod2usage(-verbose => 2, -exitval => 0) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

sub defaults {
  my ($self) = @_;
  my $o = $self->{opts};

  # Master database location:
  $o->{mhost} = 'ens-staging1' unless $o->{mhost};
  $o->{mport} = 3306 unless $o->{mport};
  
  # User database location (default values):
  $o->{port} = 3306 unless $o->{port};
  
  foreach my $required (qw/host user muser/) {
    my $msg = "Required parameter --${required} was not given";
    pod2usage(-msg => $msg, -verbose => 1, -exitval => 1) if ! $o->{$required};
  }
  
  if($o->{database} && $o->{pattern}) {
    my $msg = q{Command line options -d/--database and --pattern are mututally exclusive};
    pod2usage(-msg => $msg, -verbose => 1, -exitval => 1);
  }
  
  if($o->{pattern}) {
    my $pattern = $o->{pattern};
    $o->{pattern} = qr/$pattern/;
  }
  
  return;
}

sub databases {
  my ($self) = @_;
  my $dbc = $self->_source_dbc();
  my $o = $self->{opts};
  my ($sql, $params) = ('show databases', []);
  my $database = $o->{database};
  if($database) {
    $self->v('Filtering databases with the SQL like expression %s', $database);
    $sql .= ' like ?';
    $params = [$database];
  }
  my $dbs = $dbc->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => $params);
  my $pattern = $o->{pattern};
  if($pattern) {
    $self->v('Filtering databases with the pattern %s', $pattern);
    $dbs = [ grep { $_ =~ $pattern } @{$dbs} ];
  }
  $self->v('Found %d database(s)', scalar(@{$dbs}));
  return $dbs;
}

sub process_db {
  my ($self, $dbname) = @_;
  $self->v('Processing %s', $dbname);
  $self->v('Querying species meta');
  my $source_meta = $self->_get_meta($dbname);
  $self->v('Querying production');
  my $production = $self->_get_production($source_meta);
  $self->v('Updating attributes');
  $self->_update_attributes($source_meta, $production);
  $self->v('Updating species aliases');
  $self->_update_species_aliases($source_meta, $production);
  $self->v('Done');
  return;
}

sub _get_meta {
  my ($self, $dbname) = @_;
  my $source_dba = $self->_source_dba($dbname);
  my $mc = $source_dba->get_MetaContainer();
  my %output;
  my %tmp_alias_hash = map {lc($_) => 1} @{$mc->list_value_by_key('species.alias')};
  $output{active_aliases} = [keys %tmp_alias_hash];
  $output{inactive_aliases} = [];
  $output{production_name} = $mc->get_production_name();
  $output{scientific_name} = $mc->get_scientific_name();
  
  #Find url_name or make it from production_name
  my $url_name = $mc->list_value_by_key('species.url');
  if(@{$url_name}) {
    $output{url_name} = $url_name->[0];
  }
  else {
    $output{url_name} = ucfirst($output{production_name});
  }
  
  #Find display_name or use ensembl_alias_name ... goes into web_name in production
  my $display_name = $mc->list_value_by_key('species.display_name');
  if(@{$display_name}) {
    $output{web_name} = $display_name->[0];
  }
  else {
    $output{web_name} = $mc->single_value_by_key('species.ensembl_alias_name'); # will be replaced by just display_name
  }
  
  #Extract DB name
  my ($db_name) = $dbname =~ /^([a-z0-9_]+)_core/;
  $output{db_name} = $db_name;
  
  #Common name
  $output{common_name} = $mc->single_value_by_key('species.common_name');
  
  #Taxon
  $output{taxon} = $mc->get_taxonomy_id();
  
  return \%output;
}

sub _get_production {
  my ($self, $source_hash) = @_;
  my $target_dbc = $self->_target_dbc();
  my $hash;
  $target_dbc->sql_helper()->execute_no_return(
    -SQL => 'select species_id, common_name, web_name, scientific_name, production_name, url_name, taxon from species where taxon =?',
    -PARAMS => [$source_hash->{taxon}],
    -CALLBACK => sub {
      my ($r) = @_;
      my $i = 0;
      $hash = { 
        species_id => $r->[$i++], 
        common_name => $r->[$i++], 
        web_name => $r->[$i++], 
        scientific_name => $r->[$i++], 
        production_name => $r->[$i++], 
        url_name => $r->[$i++], 
        taxon => $r->[$i++]
      };
    }
  );
  my $sql = 'select alias from species_alias where species_id =? and is_current = ?';
  $hash->{active_aliases} = $target_dbc->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => [$hash->{species_id}, 1]);
  $hash->{inactive_aliases} = $target_dbc->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => [$hash->{species_id}, 0]);
  return $hash;
}

sub _decipher_aliases {
  my ($self, $source_meta, $production) = @_;
  
  my %source_active_hash = map { $_ => 1 } @{$source_meta->{active_aliases}};
  my %prod_active_hash = map { $_ => 1 } @{$production->{active_aliases}};
  my %prod_inactive_hash = map { $_ => 1 } @{$production->{inactive_aliases}};
  
  my @new_aliases;
  my @reactivated_aliases;
  my @retired_aliases;
  foreach my $source (keys %source_active_hash) {
    if(! $prod_active_hash{$source}) {
      if(! $prod_inactive_hash{$source}) {
        push(@new_aliases, $source);
      }
      else {
        push(@reactivated_aliases, $source);
      }
    }
  }
  foreach my $target (keys %prod_active_hash) {
    if(! $source_active_hash{$target}) {
      push(@retired_aliases, $target);
    }
  }
  
  return { new => \@new_aliases, retired => \@retired_aliases, reactivated => \@reactivated_aliases };
}

sub _decipher_attributes {
  my ($self, $source, $target) = @_;
  my %final;
  my @keys = qw/common_name web_name scientific_name production_name url_name taxon/;
  foreach my $key (@keys) {
    if($target->{$key}) {
      $final{$key} = $target->{$key};
    }
    else {
      $final{$key} = $source->{$key};
    }
  }
  return \%final;
}

sub _update_attributes {
  my ($self, $source_meta, $production) = @_;
  my $target_dbc = $self->_target_dbc();
  
  my $attributes = $self->_decipher_attributes($source_meta, $production);
    
  my @keys = keys %{$attributes};
  my $set_string = join(q{,}, map { "$_ =?" } @keys );
  my @params = map { $attributes->{$_} } @keys;
  push(@params, $production->{species_id});
  my $sql = sprintf('update species set %s where species_id =?', $set_string);
  $target_dbc->sql_helper()->execute_update(-SQL => $sql, -PARAMS => \@params);
  
  return;
}

sub _update_species_aliases {
  my ($self, $source_meta, $production) = @_;
  my $target_dbc = $self->_target_dbc();
  my $h = $target_dbc->sql_helper();

  my $aliases = $self->_decipher_aliases($source_meta, $production);
  my $species_id = $production->{species_id};
  
  foreach my $alias (@{$aliases->{new}}) {
    my $sql = 'insert into species_alias (species_id, alias, is_current, created_by, created_at) values (?,?,?,?,now())';
    $h->execute_update(-SQL => $sql, -PARAMS => [$species_id, $alias, 1, undef]);
  }
  $self->v('Added %d aliases', scalar(@{$aliases->{new}}));
  
  foreach my $alias (@{$aliases->{reactivated}}) {
    my $sql = 'update species_alias set is_current =?, modified_by =?, modified_at =now() where species_id =? and alias =?';
    $h->execute_update(-SQL => $sql, -PARAMS => [1, undef, $species_id, $alias]);
  }
  $self->v('Reactivated %d aliases', scalar(@{$aliases->{reactivated}}));
  
  if($self->{opts}->{retirealiases}){
    foreach my $alias (@{$aliases->{retired}}) {
      my $sql = 'update species_alias set is_current =?, modified_by =?, modified_at =now() where species_id =? and alias =?';
      $h->execute_update(-SQL => $sql, -PARAMS => [0, undef, $species_id, $alias]);
    }
    $self->v('Retired %d aliases', scalar(@{$aliases->{retired}}));
  }
  
  return;  
}

sub _source_dba {
  my ($self, $dbname) = @_;
  my $o = $self->{opts};
  my $dbc = $self->_source_dbc($dbname);
  return Bio::EnsEMBL::DBSQL::DBAdaptor->new(-DBCONN => $dbc, -SPECIES => $dbname);
}

sub _source_dbc {
  my ($self, $dbname) = @_;
  my $o = $self->{opts};
  my %args = ( -HOST => $o->{host}, -PORT => $o->{port}, -USER => $o->{user} );
  $args{-PASS} = $o->{pass} if $o->{pass};
  $args{-DBNAME} = $dbname if $dbname;
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub _target_dbc {
  my ($self) = @_;
  my $o = $self->{opts};
  my %args = ( -HOST => $o->{mhost}, -PORT => $o->{mport}, -DBNAME => $o->{mdatabase}, -USER => $o->{muser} );
  $args{-PASS} = $o->{mpass} if $o->{mpass};
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub v {
  my ($self, $msg, @args) = @_;
  return unless $self->{opts}->{verbose};
  my $s_msg = sprintf($msg, @args);
  my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) =
    localtime(time());
  print sprintf("[%02d-%02d-%04d %02d:%02d:%02d] %s\n",
                $mday, $mon, $year + 1900,
                $hour, $min, $sec, $s_msg);
  return;
}

__PACKAGE__->run();

__END__

=pod 

=head1 NAME

import_species_names_from_cores.pl

=head1 SYNOPSIS

  ./import_species_names_from_cores.pl -h host [-P port] \\
    -u user [-p password]
    -d database | --pattern pattern \\
    [-mh host] [-mP port] \\
    [-mu user] [-mp password] [-md database] \\
    [-v]

=head1 DESCRIPTION

This script brings in the current values for all additional fields
in the production species table such as C<production_name>, C<web_name> etc. 
It should be used as a way of synchronising these changes but in reality
these changes should be performed at the admin interface level.

=head1 OPTIONS

=over 8

=item B<-h|--host>

Host of the core databases

=item B<-P|--port>

Port of the core databases

=item B<-u|--user>

Username for the core databases

=item B<-p|--pass>

Password for the core databases

=item B<-mh|--mhost>

Host for the production database

=item B<-mP|--mport>

Port for the production database

=item B<-mu|--muser>

User for the production database

=item B<-mp|--mpass>

Pass for the production database

=item B<-md|--mdatabase>

Database name for the production database

=item B<--pattern>

Pattern to apply to filter databases by. Should be a regular expression

  --pattern="^homo.*(rnaseq|vega)_62"

=item B<-d|--database>

Database name to search for. Can be a SQL like statement

  --database="homo_sapiens_core_65_37"
  --database="%core_65%"

=item B<--retirealiases>

If specified this will retire aliases if they were not an active name in the
core database meta aliases.

=item B<--verbose>

Make the script chatty

=item B<--help>

Help message

=item B<--man>

Man page

=back

=cut
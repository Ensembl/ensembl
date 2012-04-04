#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw( :config no_ignore_case );
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Scalar qw/wrap_array/;

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->defaults();
  my $databases = $self->databases();
  foreach my $db (@{$databases}) {
    $self->v('Processing %s', $db);
    $self->_meta($db);  
    $self->v('Done');
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
    thost|th=s
    tport|tP=i
    tuser|tu=s
    tpass|tp=s
    tdatabase|td=s
    host|h=s
    port|P=i
    user|u=s
    pass|p=s
    database|d=s
    pattern=s
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
  $o->{muser} = 'ensro' unless $o->{muser};
  $o->{mdatabase} = 'ensembl_production' unless $o->{mdatabase};
  
  # Taxonomy database location:
  $o->{thost} = 'ens-staging1' unless $o->{thost};
  $o->{tuser} = 'ensro' unless $o->{tuser};
  $o->{tport} = 3306 unless $o->{tport};
  
  # User database location (default values):
  $o->{port} = 3306 unless $o->{port};
  
  foreach my $required (qw/host user tdatabase/) {
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

sub _production {
  my ($self, $db) = @_;
  $self->v('Querying production');
  my $h = $self->_production_dbc()->sql_helper();
  my $taxon = $self->_db_to_taxon($db);
  my $hash = { 'species.taxonomy_id' => $taxon };
  my $sql = 'select common_name, web_name, scientific_name, production_name, url_name, species_prefix from species where taxon =?';
  $h->execute_no_return(-SQL => $sql, -PARAMS => [$taxon], -CALLBACK => sub {
    my ($row) = @_;
    my $i = 0;
    $hash->{'species.common_name'} = $row->[$i++];
    $hash->{'species.display_name'} = $row->[$i++];
    $hash->{'species.scientific_name'} = $row->[$i++];
    $hash->{'species.production_name'} = $row->[$i++];
    $hash->{'species.url'} = $row->[$i++];
    $hash->{'species.stable_id_prefix'} = $row->[$i++];
    return;
  });
  $hash->{'species.alias'} = $h->execute_simple(
    -SQL => 'select sa.alias from species_alias sa join species s using (species_id) where s.taxon =?',
    -PARAMS => [$taxon]
  );
  return $hash;
}

sub _taxonomy {
  my ($self, $db) = @_;
  $self->v('Querying taxonomy for classification');
  my $taxon = $self->_db_to_taxon($db);
  my $excluded_ranks = [qw/species root genus/];
  my $sql = <<'SQL';
select n.name
from ncbi_taxa_node t 
join ncbi_taxa_node t2 on (t2.left_index <= t.left_index and t2.right_index >= t.right_index)
join ncbi_taxa_name n on (t2.taxon_id = n.taxon_id)
where t.taxon_id =?
and n.name_class =? 
and t2.tank not in (?,?,?)
order by t2.left_index desc
SQL
  my $res = $self->_taxon_dbc()->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => [$taxon, 'scientific name', $excluded_ranks]);
  $self->v('Classification is [%s]', join(q{, }, @{$res}));
  return $res;
}

sub _meta {
  my ($self, $db) = @_;
  my $production = $self->_production($db);
  my $taxonomy = $self->_taxonomy($db);
  $production->{'species.taxonomy_id'} = $taxonomy;
  
  $self->v('Updating meta');
  my $dba = $self->_core_dba($db);
  my $mc = $dba->get_MetaContainer();
  foreach my $key (keys %{$production}) {
    my $array = wrap_array($production->{$key});
    $mc->delete_key($key);
    foreach my $value (@{$array}) {
      $mc->store_key_value($key, $value);
    }
  }
  return;
}

sub _db_to_taxon {
  my ($self, $db) = @_;
  
  #Look at cache
  my $taxon = $self->{_db_to_taxon}->{$db};
  return $taxon if $taxon;
  
  #Try DB first
  $taxon = $self->_core_dba($db)->get_MetaContainer()->single_value_by_key('species.taxonomy_id');
  if(!$taxon) {
    die "Cannot discover the taxonomy id for the database $db. Populate meta with 'species.taxonomy_id'";
  }
  $self->{_db_to_taxon}->{$db} = $taxon;
  return $taxon;
}

sub _core_dba {
  my ($self, $dbname) = @_;
  my $o = $self->{opts};
  my $dbc = $self->_core_dbc($dbname);
  return Bio::EnsEMBL::DBSQL::DBAdaptor->new(-DBCONN => $dbc, -SPECIES => $dbname);
}

sub _core_dbc {
  my ($self, $dbname) = @_;
  my $o = $self->{opts};
  my %args = ( -HOST => $o->{host}, -PORT => $o->{port}, -USER => $o->{user} );
  $args{-PASS} = $o->{pass} if $o->{pass};
  $args{-DBNAME} = $dbname if $dbname;
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub _production_dbc {
  my ($self) = @_;
  my $o = $self->{opts};
  my %args = ( -HOST => $o->{mhost}, -PORT => $o->{mport}, -DBNAME => $o->{mdatabase}, -USER => $o->{muser} );
  $args{-PASS} = $o->{mpass} if $o->{mpass};
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub _taxon_dbc {
  my ($self) = @_;
  my $o = $self->{opts};
  my %args = ( -HOST => $o->{thost}, -PORT => $o->{tport}, -DBNAME => $o->{tdatabase}, -USER => $o->{tuser} );
  $args{-PASS} = $o->{tpass} if $o->{tpass};
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

populate_species_meta.pl

=head1 SYNOPSIS

  ./populate_species_meta.pl -h host [-P port] \\
    -u user [-p password]
    -d database \\
    [-mh host] [-mP port] \\
    [-mu user] [-mp password] [-md database] \\
    [-th host] [-tP port] \\
    [-tu user] [-tp password] [-td database] \\
    [-v]

=head1 DESCRIPTION

TODO

=head1 OPTIONS

=over 8

=item B<-h|--host>

Host of the core database(s)

=item B<-P|--port>

Port of the core database(s)

=item B<-u|--user>

Username for the core database(s)

=item B<-p|--pass>

Password for the core database(s)

=item B<-mh|--mhost>

Host for the production database

=item B<-mP|--mport>

Port for the production database

=item B<-mu|--muser>

User for the production database

=item B<-mp|--mpass>

Pass for the production database

=item B<-td|--tdatabase>

Database name for the production database

=item B<-th|--thost>

Host for the taxonomy database

=item B<-tP|--tport>

Port for the taxonomy database

=item B<-tu|--tuser>

User for the taxonomy database

=item B<-tp|--tpass>

Pass for the taxonomy database

=item B<-tp|--tpass>

Pass for the taxonomy database

=item B<-td|--tdatabase>

Database name for the taxonomy database

=item B<--pattern>

Pattern to apply to filter databases by. Should be a regular expression

  --pattern="^homo.*(rnaseq|vega)_62"

=item B<-d|--database>

Database name to search for. Can be a SQL like statement

  --database="homo_sapiens_core_65_37"
  --database="%core_65%"

=item B<--verbose>

Make the script chatty

=item B<--help>

Help message

=item B<--man>

Man page

=back

=cut
#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw( :config no_ignore_case );
use Pod::Usage;
use POSIX;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Scalar qw/wrap_array/;

sub run {
  my ($class) = @_;
  my $self = bless( {}, $class );
  $self->args();
  $self->check_opts();
  my $databases = $self->databases();
  foreach my $db ( @{$databases} ) {
    $self->v( 'Processing %s', $db );
    my $backup_table = $self->_backup($db);
    $self->_meta($db);
    $self->_remove_deprecated($db);
    $self->_remove_backup($backup_table, $db);
    $self->v('Done');
  }
  return;
}

sub args {
  my ($self) = @_;

  my $opts = {

    # Master database location:
    mhost     => 'ens-staging1',
    mport     => 3306,
    muser     => 'ensro',
    mdatabase => 'ensembl_production',

    # Taxonomy database location:
    thost => 'ens-livemirror',
    tuser => 'ensro',
    tport => 3306,
    tdatabase => 'ncbi_taxonomy',

    # User database location (default values):
    port => 3306,
    
    removedeprecated => 0, 
  };

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
    dropbaks|dB!
    dumppath|dp=s
    removedeprecated!
    help
    man
    /;
  GetOptions( $opts, @cmd_opts ) or pod2usage( -verbose => 1, -exitval => 1 );
  pod2usage( -verbose => 1, -exitval => 0 ) if $opts->{help};
  pod2usage( -verbose => 2, -exitval => 0 ) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

sub check_opts {
  my ($self) = @_;
  my $o = $self->{opts};

  foreach my $required (qw/host user tdatabase/) {
    my $msg = "Required parameter --${required} was not given";
    pod2usage( -msg => $msg, -verbose => 1, -exitval => 1 ) if !$o->{$required};
  }

  if ( $o->{database} && $o->{pattern} ) {
    my $msg =
q{Command line options -d/--database and --pattern are mututally exclusive};
    pod2usage( -msg => $msg, -verbose => 1, -exitval => 1 );
  }

  if ( $o->{pattern} ) {
    my $pattern = $o->{pattern};
    $o->{pattern} = qr/$pattern/;
  }
  
  if($o->{dropbaks} && ! defined $o->{dumppath} ) {
    die "If you are specifying --dropbaks you must specify a --dumppath";
  }
  if(defined $o->{dumppath} && ! -d $o->{dumppath}) {
    die sprintf('--dumppath %s does not exist or is not a directory', $o->{dumppath});
  }

  return;
}

sub databases {
  my ($self) = @_;
  my $dbc    = $self->_core_dbc();
  my $o      = $self->{opts};
  my ( $sql, $params ) = ( 'show databases', [] );
  my $database = $o->{database};
  if ($database) {
    $self->v( 'Filtering databases with the SQL like expression %s',
      $database );
    $sql .= ' like ?';
    $params = [$database];
  }
  my $dbs =
    $dbc->sql_helper()->execute_simple( -SQL => $sql, -PARAMS => $params );
  my $pattern = $o->{pattern};
  if ($pattern) {
    $self->v( 'Filtering databases with the pattern %s', $pattern );
    $dbs = [ grep { $_ =~ $pattern } @{$dbs} ];
  }
  $self->v( 'Found %d database(s)', scalar( @{$dbs} ) );
  return $dbs;
}

sub _production {
  my ( $self, $db ) = @_;
  $self->v('Querying production');
  my $dbc   = $self->_production_dbc();
  my $h     = $dbc->sql_helper();
  my $taxon = $self->_db_to_taxon($db);
  my $hash  = { 'species.taxonomy_id' => $taxon };
  my $sql =
'select common_name, web_name, scientific_name, production_name, url_name, species_prefix from species where taxon =?';
  $h->execute_no_return(
    -SQL      => $sql,
    -PARAMS   => [$taxon],
    -CALLBACK => sub {
      my ($row) = @_;
      my $i = 0;
      $hash->{'species.common_name'}      = $row->[ $i++ ];
      $hash->{'species.display_name'}     = $row->[ $i++ ];
      $hash->{'species.scientific_name'}  = $row->[ $i++ ];
      $hash->{'species.production_name'}  = $row->[ $i++ ];
      $hash->{'species.url'}              = $row->[ $i++ ];
      $hash->{'species.stable_id_prefix'} = $row->[ $i++ ];
      return;
    }
  );
  $hash->{'species.alias'} = $h->execute_simple(
    -SQL =>
'select sa.alias from species_alias sa join species s using (species_id) where s.taxon =?',
    -PARAMS => [$taxon]
  );
  return $hash;
}

sub _taxonomy {
  my ( $self, $db ) = @_;
  $self->v('Querying taxonomy for classification');
  my $taxon          = $self->_db_to_taxon($db);
  my @excluded_ranks = ('root', 'genus', 'species subgroup', 'species group', 'subgenus');
  my @excluded_names = ('cellular organisms', 'root');
  my $sql            = <<'SQL';
select n.name
from ncbi_taxa_node t 
join ncbi_taxa_node t2 on (t2.left_index <= t.left_index and t2.right_index >= t.right_index)
join ncbi_taxa_name n on (t2.taxon_id = n.taxon_id)
where t.taxon_id =?
and n.name_class =? 
and t2.rank not in (?,?,?,?,?)
and n.name not in (?,?)
order by t2.left_index desc
SQL
  my $dbc = $self->_taxon_dbc();
  my $res = $dbc->sql_helper()->execute_simple(
    -SQL    => $sql,
    -PARAMS => [ $taxon, 'scientific name', @excluded_ranks, @excluded_names ]
  );
  $self->v( 'Classification is [%s]', join( q{, }, @{$res} ) );
  return $res;
}

sub _meta {
  my ( $self, $db ) = @_;
  my $production = $self->_production($db);
  my $taxonomy   = $self->_taxonomy($db);
  $production->{'species.classification'} = $taxonomy;

  $self->v('Updating meta');
  my $dba = $self->_core_dba($db);
  my $mc  = $dba->get_MetaContainer();
  foreach my $key ( keys %{$production} ) {
    my $array = wrap_array( $production->{$key} );
    $mc->delete_key($key);
    foreach my $value ( @{$array} ) {
      $mc->store_key_value( $key, $value );
    }
  }
  return;
}

sub _backup {
  my ( $self, $dbname ) = @_;
  my $o = $self->{opts};

  my $core = $self->_core_dbc($dbname);
  my $table = 'meta_bak';
  $core->do('drop table if exists '.$table);
  
  $self->v( 'Backing up to %s', $table );
  $core->do( sprintf( 'create table %s like meta', $table ) );
  $self->v( 'Copying data from meta to %s', $table );
  $core->do( sprintf( 'insert into %s select * from meta', $table ) );
  $self->v('Done backup');

  if ( defined($o->{dumppath}) ) {
    my $timestamp = strftime( "%Y%m%d-%H%M%S", localtime() );
    
    # Backup the table on file.
    my $filename = sprintf( "%s/%s.%s.%s.sql",
                            $o->{dumppath}, $dbname,
                            $table,    $timestamp );

    if ( -e $filename ) {
      die( sprintf( "File '%s' already exists.", $filename ) );
    }

    $self->v( "Backing up table %s on file", $table );
    $self->v( "--> %s",                       $filename );

    if (system( join(q{ }, "mysqldump",
                "--host=".$o->{host},
                "--port=".$o->{port},
                "--user=".$o->{user},
                (
                  defined($o->{pass}) ? "--password=".$o->{pass} : ""
                ),
                "--result-file=$filename",
                $dbname,
                $table)))
    {
      die("mysqldump failed: $?");
    }
  }

  return $table;
}

sub _remove_deprecated {
  my ($self, $db) = @_;
  
  if(!$self->{opts}->{removedeprecated}) {
    $self->v('Not removing deprecated meta keys');
    return;
  }
  
  my $dba = $self->_core_dba($db);
  my $mc  = $dba->get_MetaContainer();
  
  $self->v('Removing deprecated meta keys');
  my @deprecated_keys = qw/
    species.ensembl_common_name species.ensembl_alias_name
    species.short_name
  /;
  foreach my $key (@deprecated_keys) {
    $self->v('Deleting key "%s"', $key);
    $mc->delete_key($key);
  }
  $self->v('Finished removing deprecated meta keys');
  
  return;
}

sub _remove_backup {
  my ($self, $table, $dbname) = @_;
  if($self->{opts}->{dropbaks}) {
    $self->v('Dropping backup table %s', $table);
    my $dbc = $self->_core_dbc($dbname);
    $dbc->do('drop table if exists '.$table);
  }
  else {
    $self->v('Leaving backup table %s', $table);
  }
  return;
}

sub _db_to_taxon {
  my ( $self, $db ) = @_;

  #Look at cache
  my $taxon = $self->{_db_to_taxon}->{$db};
  return $taxon if $taxon;

  #Try DB first
  $taxon =
    $self->_core_dba($db)->get_MetaContainer()
    ->single_value_by_key('species.taxonomy_id');
  if ( !$taxon ) {
    die
"Cannot discover the taxonomy id for the database $db. Populate meta with 'species.taxonomy_id'";
  }
  $self->{_db_to_taxon}->{$db} = $taxon;
  return $taxon;
}

sub _core_dba {
  my ( $self, $dbname ) = @_;
  my $o   = $self->{opts};
  my $dbc = $self->_core_dbc($dbname);
  return Bio::EnsEMBL::DBSQL::DBAdaptor->new( -DBCONN => $dbc,
    -SPECIES => $dbname );
}

sub _core_dbc {
  my ( $self, $dbname ) = @_;
  my $o = $self->{opts};
  my %args = ( -HOST => $o->{host}, -PORT => $o->{port}, -USER => $o->{user} );
  $args{-PASS}   = $o->{pass} if $o->{pass};
  $args{-DBNAME} = $dbname    if $dbname;
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub _production_dbc {
  my ($self) = @_;
  my $o      = $self->{opts};
  my %args   = (
    -HOST   => $o->{mhost},
    -PORT   => $o->{mport},
    -DBNAME => $o->{mdatabase},
    -USER   => $o->{muser}
  );
  $args{-PASS} = $o->{mpass} if $o->{mpass};
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub _taxon_dbc {
  my ($self) = @_;
  my $o      = $self->{opts};
  my %args   = (
    -HOST   => $o->{thost},
    -PORT   => $o->{tport},
    -DBNAME => $o->{tdatabase},
    -USER   => $o->{tuser}
  );
  $args{-PASS} = $o->{tpass} if $o->{tpass};
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub v {
  my ( $self, $msg, @args ) = @_;
  return unless $self->{opts}->{verbose};
  my $s_msg = sprintf( $msg, @args );
  my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
    localtime( time() );
  print sprintf(
    "[%02d-%02d-%04d %02d:%02d:%02d] %s\n",
    $mday, $mon, $year + 1900,
    $hour, $min, $sec, $s_msg
  );
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
    [-dropbaks]
    [-dumppath]
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

=item B<-th|--thost>

Host for the taxonomy database

=item B<-tP|--tport>

Port for the taxonomy database

=item B<-tu|--tuser>

User for the taxonomy database

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

=item B<--dropbaks>

Forces the dropping of the _bak tables created as part
of this script's run. Normally we would suggest you keep
this option off unless you know you do not want these. You
must run this option with --dumppath to avoid unintentional
data loss.

=item B<--dumppath>

Back-up tables into the specified directory path. If you are
using --dropbaks you must specify this option

=item B<--removedeprecated>

Deletes any key from the DB which is deemed to be deprecated. A temporary
option until we finish working with the latest set of meta key changes.

=item B<--verbose>

Make the script chatty

=item B<--help>

Help message

=item B<--man>

Man page

=back

=cut

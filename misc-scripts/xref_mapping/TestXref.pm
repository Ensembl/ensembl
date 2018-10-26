package TestXref;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::IO qw/slurp/;
use DBI;

sub new {
  my ($class,$root_path,$do_not_clean_up) = @_;
  my $self = bless {}, $class;
  my $test_file = 'testdb.sq';
  if (defined $root_path) {
    $test_file = $root_path . $test_file;
  }
  $self->{db_path} = $test_file;
  $self->{dbh} = DBI->connect("dbi:SQLite:dbname=$test_file",undef,undef, { AutoCommit => 1});
  $self->{no_clean_up} = $do_not_clean_up;
  $self->_load_schema;
  return $self;
}

sub dbi {
  my $self = shift;
  return $self->{dbh};
}

sub check_direct_xref {
  my ($self,$xref) = @_;
  my $sth = $self->{dbh}->prepare('SELECT * FROM xref WHERE db_primaryacc = ?');
  $sth->execute($xref->{db_primaryacc});

}


sub _load_schema {
  my $self = shift;
  my $sql_file = '/Users/ktaylor/ensembl/ensembl/misc-scripts/xref_mapping/sql/table.sql';
  my $var;
  my $sql = slurp($sql_file);
  $sql =~ s/#.*?\n//mg;
  $sql =~ s/\/\*\*.*?\*\///mg;
  $sql =~ s/^--.*?\n//mg;
  $sql =~ s/^\s*$//mg;
  $sql =~ s/AUTO_INCREMENT//gi;
  print $sql."\n";
  $self->dbi->do($sql);
}


sub DESTROY {
  my $self = shift;
  if (! $self->{no_clean_up}) {
    unlink $self->{db_path};
  }
}

1;
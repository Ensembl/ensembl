
### Bio::EnsEMBL::DBSQL::SQL

package Bio::EnsEMBL::DBSQL::SQL;

use strict;
use DBI;
use Bio::EnsEMBL::DBSQL::SQL::mysql;
use Bio::EnsEMBL::DBSQL::SQL::oracle;
use Bio::EnsEMBL::DBSQL::SQL::sybase;
use Bio::EnsEMBL::DBSQL::SQL::StatementHandle;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use vars '@ISA';

@ISA = qw{ DBI::db };

sub new {
  my( $pkg, $dsn, $user, $password ) = @_;
    
  my $dbh = DBI->connect($dsn, $user, $password, 
                         {RaiseError => 0, PrintError => 0});

  return bless($dbh, $pkg) if($dbh);

  throw("Can't connect to SQL database with:\n"
        . "       dsn = '$dsn'\n"
        . "      user = '$user'\n"
        . "  password = '$password'\n");
}

sub prepare {
  my( $dbh, @args ) = @_;
    
  my $sth = $dbh->SUPER::prepare(@args);
  if ($sth) {
    bless($sth, 'Bio::EnsEMBL::DBSQL::SQL::StatementHandle');
    return $sth;
  }

  throw("prepare failed: '$DBI::errstr'");
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::DBSQL::SQL

=head1 DESCRIPTION

Ensembl's SQL compatability layer.

This module inherits from B<DBI::db>, overriding
the B<prepare> method to bless the created
statement handles into the
B<Bio::EnsEMBL::DBSQL::SQL::StatementHandle>
class.

Database handles created in
B<Bio::Ensembl::DBSQL::Obj> are blessed into one
of the database driver specific classes (eg:
B<Bio::EnsEMBL::DBSQL::SQL::mysql>) which inherit
from this class.

=head1 MEHTODS

=over 4

=item new

    my $dbh = Bio::EnsEMBL::DBSQL::SQL::<driver>->new($dsn, $user, $password);

Called by B<Bio::EnsEMBL::DBSQL::Obj::new>, it
creates a B<DBI::db> database handle, blesses it
into the SQL compatability driver package, and
returns it.

=item prepare

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk


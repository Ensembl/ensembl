
### Bio::EnsEMBL::DBOLD::SQL

package Bio::EnsEMBL::DBOLD::SQL;

use strict;
use DBI;
use Bio::EnsEMBL::DBOLD::SQL::mysql;
use Bio::EnsEMBL::DBOLD::SQL::oracle;
use Bio::EnsEMBL::DBOLD::SQL::sybase;
use Bio::EnsEMBL::DBOLD::SQL::StatementHandle;
use Bio::EnsEMBL::Root;
use vars '@ISA';

@ISA = qw{ DBI::db Bio::EnsEMBL::Root };

sub new {
    my( $pkg, $dsn, $user, $password ) = @_;
    
    my $dbh = DBI->connect($dsn, $user, $password, {RaiseError => 0, PrintError => 0});
    if ($dbh) {
        return bless($dbh, $pkg);
    } else {
        # Create a BioPerl object in order to throw an exception
        my $self = bless {}, $pkg;
        $self->throw("Can't connect to SQL database with:\n"
            . "       dsn = '$dsn'\n"
            . "      user = '$user'\n"
            . "  password = '$password'\n"
            );
    }
}

sub prepare {
    my( $dbh, @args ) = @_;
    
    my $sth = $dbh->SUPER::prepare(@args);
    if ($sth) {
        bless($sth, 'Bio::EnsEMBL::DBOLD::SQL::StatementHandle');
        return $sth;
    } else {
        $dbh->throw("prepare failed: '$DBI::errstr'");
    }
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::DBOLD::SQL

=head1 DESCRIPTION

Ensembl's SQL compatability layer.

This module inherits from B<DBI::db>, overriding
the B<prepare> method to bless the created
statement handles into the
B<Bio::EnsEMBL::DBOLD::SQL::StatementHandle>
class.

Database handles created in
B<Bio::Ensembl::DBOLD::Obj> are blessed into one
of the database driver specific classes (eg:
B<Bio::EnsEMBL::DBOLD::SQL::mysql>) which inherit
from this class.

=head1 MEHTODS

=over 4

=item new

    my $dbh = Bio::EnsEMBL::DBOLD::SQL::<driver>->new($dsn, $user, $password);

Called by B<Bio::EnsEMBL::DBOLD::Obj::new>, it
creates a B<DBI::db> database handle, blesses it
into the SQL compatability driver package, and
returns it.

=item prepare

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk


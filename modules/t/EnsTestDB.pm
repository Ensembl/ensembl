
=pod

=head1 NAME - EnsTestDB

=head1 SYNOPSIS

    # Add test dir to lib search path
    use lib 't';
    
    use EnsTestDB;
    
    my $ens_test = EnsTestDB->new();
    
    # Load some data into the db
    $ens_test->do_sql_file("some_data.sql");
    
    # Get an EnsEMBL db object for the test db
    my $db = $ens_test->get_DBSQL_Obj;

=head1 DESCRIPTION

This is a module used just by the EnsEMBL test
suite to create a test database for a particular
test.  Creating a new object creates a database
with a name such that it should never clash with
other users testing on the same server.  The
database is destroyed when the object goes out of
scope.

The settings, such as the server host and port,
are found in the file B<EnsTestDB.conf>.  See
B<EnsTestDB.conf.example> for an example.

=head1 METHODS

=cut

package EnsTestDB;

use strict;
use Sys::Hostname 'hostname';
use Bio::EnsEMBL::DBLoader;
use DBI;
use Carp;

#Package variable for unique database name
my $counter=0;

{
    # This is a list of possible entries in the config
    # file "EnsTestDB.conf"
    my %known_field = map {$_, 1} qw(
        driver
        host
        user
        port
        password
        schema_sql
        module
        );
    
    sub new {
        my( $pkg ) = @_;

        $counter++;
        # Get config from file, or use default values
        my $self = do 'EnsTestDB.conf' || {
            'driver'        => 'mysql',
            'host'          => 'localhost',
            'user'          => 'root',
            'port'          => '3306',
            'password'      => undef,
            'schema_sql'    => ['../sql/table.sql'],
            'module'        => 'Bio::EnsEMBL::DBSQL::Obj'
            };
        foreach my $f (keys %$self) {
            confess "Unknown config field: '$f'" unless $known_field{$f};
        }
        bless $self, $pkg;
        $self->create_db;
	
        return $self;
    }
}

sub driver {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'driver'} = $value;
    }
    return $self->{'driver'} || confess "driver not set";
}

sub host {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'host'} = $value;
    }
    return $self->{'host'} || confess "host not set";
}

sub user {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'user'} = $value;
    }
    return $self->{'user'} || confess "user not set";
}

sub port {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'port'} = $value;
    }
    return $self->{'port'};
}

sub password {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->{'password'} = $value;
    }
    return $self->{'password'};
}

sub schema_sql {
    my( $self, $value ) = @_;
    
    if ($value) {
        push(@{$self->{'schema_sql'}}, $value);
    }
    return $self->{'schema_sql'} || confess "schema_sql not set";
}

sub dbname {
    my( $self ) = @_;

    $self->{'_dbname'} ||= $self->_create_db_name();
    return $self->{'_dbname'};
}

sub module {
    my ($self, $value) = @_;
    $self->{'module'} = $value if ($value);
    return $self->{'module'};
}

sub _create_db_name {
    my( $self ) = @_;

    my $host = hostname();
    my $db_name = "_test_db_${host}_$$".$counter;
    $db_name =~ s{\W}{_}g;
    return $db_name;
}

sub create_db {
    my( $self ) = @_;
    
    ### FIXME: not portable between different drivers
    my $locator = 'dbi:'. $self->driver .':host='. $self->host .';database=mysql';
    my $db = DBI->connect(
        $locator, $self->user, $self->password, {RaiseError => 1}
        ) or confess "Can't connect to server";
    my $db_name = $self->dbname;
    $db->do("CREATE DATABASE $db_name");
    $db->disconnect;
    
    $self->do_sql_file(@{$self->schema_sql});
}

sub db_handle {
    my( $self ) = @_;
    
    unless ($self->{'_db_handle'}) {
        $self->{'_db_handle'} = DBI->connect(
            $self->test_locator, $self->user, $self->password, {RaiseError => 1}
            ) or confess "Can't connect to server";
    }
    return $self->{'_db_handle'};
}

sub test_locator {
    my( $self ) = @_;
    
    my $locator = 'dbi:'. $self->driver .':database='. $self->dbname;
    foreach my $meth (qw{ host port }) {
        if (my $value = $self->$meth()) {
            $locator .= ";$meth=$value";
        }
    }
    return $locator;
}

sub ensembl_locator {
    my( $self) = @_;
    
    my $module = ($self->module() || 'Bio::EnsEMBL::DBSQL::Obj');
    my $locator = '';
    foreach my $meth (qw{ host port dbname user }) {
        my $value = $self->$meth();
	if( !defined $value ) { next; }
        $locator .= ';' if $locator;
        $locator .= "$meth=$value";
    }
    return "$module/$locator";
}

sub get_DBSQL_Obj {
    my( $self ) = @_;
    
    my $locator = $self->ensembl_locator();
    return Bio::EnsEMBL::DBLoader->new($locator);
}

sub do_sql_file {
    my( $self, @files ) = @_;
    local *SQL;
    my $i = 0;
    my $dbh = $self->db_handle;
    
    foreach my $file (@files)
    {
        my $sql = '';
        open SQL, $file or die "Can't read SQL file '$file' : $!";
        while (<SQL>) {
            s/(#|--).*//;       # Remove comments
            next unless /\S/;   # Skip lines which are all space
            $sql .= $_;
            $sql .= ' ';
        }
        close SQL;
        
        foreach my $s (grep /\S/, split /;/, $sql) {
            $self->validate_sql($s);
            $dbh->do($s);
            $i++
        }
    }
    return $i;
}

sub validate_sql {
    my ($self, $statement) = @_;
    if ($statement =~ /insert/i)
    {
        $statement =~ s/\n/ /g; #remove newlines
        die ("INSERT should use explicit column names (-c switch in mysqldump)\n$statement\n")
            unless ($statement =~ /insert.+into.*\(.+\).+values.*\(.+\)/i);
    }
}

sub DESTROY {
    my( $self, $file ) = @_;
    
    if (my $dbh = $self->db_handle) {
        my $db_name = $self->dbname;
        $dbh->do("DROP DATABASE $db_name");
        $dbh->disconnect;
    }
}

1;


__END__

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk


### Bio::EnsEMBL::DBSQL::SQL::StatementHandle

package Bio::EnsEMBL::DBSQL::SQL::StatementHandle;

use strict;
use DBI;
use Bio::EnsEMBL::Root;
use vars '@ISA';

@ISA = qw{ DBI::st Bio::EnsEMBL::Root };


sub execute {
    my( $sth, @args ) = @_;
    
    my $result = $sth->SUPER::execute(@args);
    if ($result) {
        return $result;
    } else {
        $sth->throw("execute failed : '$DBI::errstr'");
    }
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::DBSQL::SQL::StatementHandle

=head1 DESCRIPTION

This package inherits from B<DBI::st>, and
overrides its B<execute> function to provide nice
BioPerl style exceptions.

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk


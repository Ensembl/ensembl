
### Bio::EnsEMBL::DBSQL::SQL::StatementHandle

package Bio::EnsEMBL::DBSQL::SQL::StatementHandle;

use strict;
use DBI;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use vars '@ISA';

@ISA = qw{ DBI::st };


###
### Uncomment the following to turn on useful debugging info
### you will also have to alter DBConnection.pm
###

#BEGIN {
# eval {
# require Time::HiRes;
# Time::HiRes->import('time');
# };
#};

#sub execute {
#  my( $sth, @args ) = @_;
#  my $time = time;
#  my $result = $sth->SUPER::execute(@args);
#  $time = time - $time;
#  print STDERR "query time: $time\n";
#  if ($result) {
#    return $result;
#  }  
#  throw("execute failed : '$DBI::errstr'");
#}

1;

__END__

=head1 NAME - Bio::EnsEMBL::DBSQL::SQL::StatementHandle

=head1 DESCRIPTION

This package inherits from B<DBI::st>, and
overrides its B<execute> function to provide nice
BioPerl style exceptions.

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk


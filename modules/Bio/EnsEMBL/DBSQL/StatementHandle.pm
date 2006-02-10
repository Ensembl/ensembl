=head1 NAME - Bio::EnsEMBL::DBSQL::StatementHandle

=head1 SYNOPSIS

  Do not use this class directly.  It will automatically be used by
  the Bio::EnsEMBL::DBSQL::DBConnection class.

=head1 DESCRIPTION

  This class extends DBD::mysql::st so that the DESTROY method may be
  overridden.  If the DBConnection::disconnect_when_inactive flag is set
  this statement handle will cause the database connection to be closed
  when it goes out of scope and there are no other open statement handles.

=head1 CONTACT

  This module is part of the Ensembl project: www.ensembl.org

  Ensembl development mailing list: <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::StatementHandle;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(warning stack_trace_dump);

use DBD::mysql;
use DBI;

# use Time::HiRes qw(time);

@ISA = qw(DBI::st);

# As DBD::mysql::st is a tied hash can't store things in it,
# so have to have parallel hash
my %dbchash;
my %dbc_sql_hash;
sub dbc {
  my $self = shift;

  if (@_) {
    my $dbc = shift;
    if(!defined($dbc)) {
      # without delete key space would grow indefinately causing mem-leak
      delete($dbchash{$self});
    } else {
      $dbchash{$self} = $dbc;
    }
  }

  return $dbchash{$self};
}

sub sql {
  my $self = shift;

  if (@_) {
    my $sql = shift;
    if(!defined($sql)) {
      # without delete key space would grow indefinately causing mem-leak
      delete($dbc_sql_hash{$self});
    } else {
      $dbc_sql_hash{$self} = $sql;
    }
  }

  return $dbc_sql_hash{$self};
}


#
# uncomment this for printing out handy debug information (every query)
#
#sub execute {
#   my $self = shift;
#
#   my $sql = $self->sql();
#
#
#   my @chrs = split(//, $sql);
#
#   my $j = 0;
#   for(my $i =0; $i < @chrs; $i++) {
#     $chrs[$i] = $_[$j++] if($chrs[$i] eq '?' && defined($_[$j]));
#   }
#
#   my $str = join('', @chrs);
#
#   my $time = time;
##   print STDERR "\nSQL:\n$str\n\n";
#
# #  print STDERR stack_trace_dump(), "\n";
#
#   my $res = $self->SUPER::execute(@_);
#   $time = time - $time;
##   print STDERR "DONE ($time)\n";
#
#   return $res;
# }


sub DESTROY {
  my ($obj) = @_;

  my $dbc = $obj->dbc;
  $obj->dbc(undef);
  my $sql = $obj->sql;
  $obj->sql(undef);

  # rebless into DBI::st so that superclass destroy method is called
  # if it exists (it does not exist in all DBI versions)
  bless($obj, 'DBI::st');

  # The count for the number of kids is decremented only after this
  # function is complete. Disconnect if there is 1 kid (this one) remaining.
  if($dbc  && $dbc->disconnect_when_inactive() &&
     $dbc->connected && ($dbc->db_handle->{Kids} == 1))
  {
    if($dbc->disconnect_if_idle()) {
       warn("Problem disconnect $obj around sql = $sql\n");
    }
  }
}

1;

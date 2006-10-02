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

#use Time::HiRes qw(time);

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
# uncomment this for printing out handy debug information 
# (every query if you want)
#
## call   Bio::EnsEMBL::DBSQL::StatementHandle->dump(1); to start log
## call   Bio::EnsEMBL::DBSQL::StatementHandle->dump(0); to end log
## and set $dump to 0 
## leave $dump = 1 for continuous log
#my @bind_args=();
#my $dump = 0;
#my $total_time = 0;

#sub dump {
#  $dump = shift;
#  if($total_time){
#    print STDERR "###################Time taken was $total_time\n";
#  }
#  $total_time = 0;
#}

#sub bind_param{ 
#  my ($self,@args) = @_;

#  $bind_args[$args[0]-1] = $args[1];
#  $self->SUPER::bind_param(@args);
#} 

#sub execute {
#  my ($self,@args) = @_;
  
#  if(!$dump){ # skip dumping
#    return $self->SUPER::execute(@args);
#  }
#  my $sql = $self->sql();
  
  
#  my @chrs = split(//, $sql);
  
#  my $j = 0;
#  for(my $i =0; $i < @chrs; $i++) {
#    $chrs[$i] = $bind_args[$j++] if($chrs[$i] eq '?' && defined($bind_args[$j]));
#  }
  
#  my $str = join('', @chrs);
  
#  print STDERR "\nSQL:\n$str\n\n";
  
#  my $time = time;
#  my $res = $self->SUPER::execute(@args);
#  $time = time - $time;

#  $total_time += $time;
#  print "DONE ($time)\n";
#  return $res;
#}

# End uncomment


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

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
      # without delete key space would grow indefinitely causing mem-leak
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
      # without delete key space would grow indefinitely causing mem-leak
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
#
# to stop caching messing up your timings try doing thefoloowing on any adapter:-
#
# my $sth = $slice_adaptor->dbc->db_handle->prepare("SET SESSION query_cache_type = OFF");
# $sth->execute || die "set session failed\n";
# 
#
## call   Bio::EnsEMBL::DBSQL::StatementHandle->sql_timing_start(); to start log
## call   Bio::EnsEMBL::DBSQL::StatementHandle->sql_timing_print(1); to print the results
## and set $dump to 0 
## leave $dump = 1 for continuous log
# uncomment from here-------------------
#my @bind_args=();
#my $dump = 0;
#my %total_time;
#my %min_time;
#my %max_time;
#my %number_of_times;
#my %first_time;
#my $grand_total;

#sub sql_timing_start{
#  %total_time = ();
#  %number_of_times = ();
#  %min_time = ();
#  %max_time = ();
#  %first_time = ();
#  $grand_total = 0;
#  $dump = 1;
#}

#sub sql_timimg_pause{
#  $dump=0;
#}

#sub sql_timing_resume{
#  $dump =1;
#}



#sub sql_timing_print{
#  my $self = shift;
#  my $level = shift;
#  my $fh    = shift;
#  my $grand_total=0;

#  if(!defined($fh)){
#    $fh = \*STDERR;
#  }
#  print ref($fh)."\n";
#  foreach my $key (keys %total_time){
#    if(defined($level) and $level){
#      print $fh  $key."\n";
#      print $fh "total\t \tnum\tfirst \t\tavg\t \t[min     ,max      ]\n";
#      printf $fh "%6f\t%d\t%6f\t%6f\t[%6f, %6f]\n\n", 
#	$total_time{$key}, 
#	$number_of_times{$key},
#        $first_time{$key},
#       ($total_time{$key}/$number_of_times{$key}),
#        $min_time{$key},
#	$max_time{$key};
#    }	
#    $grand_total += $total_time{$key};
#  }
#  printf $fh "\ntotal time %6f\n\n", $grand_total;
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
  
##  uncomment this line if you want to see sql in order.
##  print STDERR "\n\nSQL:\n$str\n\n";
  
#  my $time = time;
#  my $res = $self->SUPER::execute(@args);
#  $time = time - $time;

#  if(defined($total_time{$sql})){
#    $total_time{$sql} += $time;
#    $number_of_times{$sql}++;
#    if($min_time{$sql} > $time){
#      $min_time{$sql} = $time;
#    }
#    if($max_time{$sql} < $time){
#      $max_time{$sql} = $time;
#    }
#  }
#  else{
#    $first_time{$sql} = $time;
#    $max_time{$sql} = $time;
#    $min_time{$sql} = $time;
#    $total_time{$sql} = $time;
#    $number_of_times{$sql} = 1;
#  }
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

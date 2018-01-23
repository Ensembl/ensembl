=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::StatementHandle

=head1 SYNOPSIS

Do not use this class directly.  It will automatically be used by the
Bio::EnsEMBL::DBSQL::DBConnection class.

=head1 DESCRIPTION

This class extends DBD::mysql::st so that the DESTROY method may be
overridden.  If the DBConnection::disconnect_when_inactive flag is set
this statement handle will cause the database connection to be closed
when it goes out of scope and there are no other open statement handles.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::StatementHandle;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(warning throw);

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

sub DESTROY {
  my ($self) = @_;

  my $dbc = $self->dbc;
  $self->dbc(undef);
  my $sql = $self->sql;
  $self->sql(undef);

  # Re-bless into DBI::st so that superclass destroy method is called if
  # it exists (it does not exist in all DBI versions).
  bless( $self, 'DBI::st' );

  # The count for the number of kids is decremented only after this
  # function is complete. Disconnect if there is 1 kid (this one)
  # remaining.
  if (    $dbc
       && $dbc->disconnect_when_inactive()
       && $dbc->connected
       && ( $dbc->db_handle->{Kids} == 1 ) )
  {
    if ( $dbc->disconnect_if_idle() ) {
      warn("Problem disconnect $self around sql = $sql\n");
    }
  }
} ## end sub DESTROY

1;

# Comment out this "__END__" for printing out handy debug information
# (every query if you want).

__END__

# To stop caching messing up your timings, try doing the following on
# any adapter:
#
#   $slice_adaptor->dbc()->db_handle()
#       ->do("SET SESSION query_cache_type = OFF");
#
# To start logging:
# Bio::EnsEMBL::DBSQL::StatementHandle->sql_timing_start();
#
# To display the results:
# Bio::EnsEMBL::DBSQL::StatementHandle->sql_timing_print(1);
#
# To pause logging:
# Bio::EnsEMBL::DBSQL::StatementHandle->sql_timimg_pause();
#
# To resume logging after pause:
# Bio::EnsEMBL::DBSQL::StatementHandle->sql_timimg_resume();

use Time::HiRes qw(time);

my @bind_args = ();
my $dump      = 0;
my %total_time;
my %min_time;
my %max_time;
my %number_of_times;
my %first_time;
my $grand_total;

sub sql_timing_start {
  %total_time      = ();
  %number_of_times = ();
  %min_time        = ();
  %max_time        = ();
  %first_time      = ();
  $dump            = 1;
}

sub sql_timing_pause  { $dump = 0 }
sub sql_timing_resume { $dump = 1 }

sub sql_timing_print {
  my ( $self, $level, $fh ) = @_;

  my $grand_total = 0;

  if ( !defined($fh) ) {
    $fh = \*STDERR;
  }

  print( ref($fh), "\n" );

  foreach my $key ( keys %total_time ) {
    $grand_total += $total_time{$key};

    if ( !( defined($level) and $level ) ) { next }

    print( $fh $key, "\n" );

    print( $fh
        "total\t \tnum\tfirst \t\tavg\t \t[min     ,max      ]\n" );

    printf( $fh "%6f\t%d\t%6f\t%6f\t[%6f, %6f]\n\n",
      $total_time{$key}, $number_of_times{$key},
      $first_time{$key}, ( $total_time{$key}/$number_of_times{$key} ),
      $min_time{$key}, $max_time{$key} );
  }

  printf( $fh "\ntotal time %6f\n\n", $grand_total );

} ## end sub sql_timing_print

sub bind_param {
  my ( $self, @args ) = @_;

  $bind_args[ $args[0] - 1 ] = $args[1];
  $self->SUPER::bind_param(@args);
}

sub execute {
  my ( $self, @args ) = @_;

  my $retval;
  # Skip dumping if !$dump
  if ( !$dump ) {
      local $self->{RaiseError};
      $retval = $self->SUPER::execute(@args);
      if ( !defined($retval) ) {
        throw("Failed to execute SQL statement");
      }
      return $retval;
  }

  my $sql = $self->sql();
  my @chrs = split( //, $sql );

  my $j = 0;

  for ( my $i = 0; $i < @chrs; $i++ ) {
    if ( $chrs[$i] eq '?' && defined( $bind_args[$j] ) ) {
      $chrs[$i] = $bind_args[ $j++ ];
    }
  }

  my $str = join( '', @chrs );

  # Uncomment this line if you want to see sql in order.
  # print( STDERR "\n\nSQL:\n$str\n\n" );

  my $time = time();
  {
    local $self->{RaiseError};
    $retval = $self->SUPER::execute(@args);
    if ( !defined($retval) ) {
      throw("Failed to execute SQL statement");
    }
  }
  #  my $res  = $self->SUPER::execute(@args);
  $time = time() - $time;

  if ( defined( $total_time{$sql} ) ) {
    $total_time{$sql} += $time;
    $number_of_times{$sql}++;

    if ( $min_time{$sql} > $time ) { $min_time{$sql} = $time }
    if ( $max_time{$sql} < $time ) { $max_time{$sql} = $time }

  } else {
    $first_time{$sql}      = $time;
    $max_time{$sql}        = $time;
    $min_time{$sql}        = $time;
    $total_time{$sql}      = $time;
    $number_of_times{$sql} = 1;
  }

  return $retval;
} ## end sub execute

1;

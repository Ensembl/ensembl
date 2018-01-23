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

Bio::EnsEMBL::DBSQL::ProxyDBConnection - Database connection wrapper allowing
for one backing connection to be used for multiple DBs

=head1 SYNOPSIS

  my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(-HOST => 'host', -PORT => 3306, -USER => 'user');
  my $p_h_dbc = Bio::EnsEMBL::DBSQL::ProxyDBConnection->new(-DBC => $dbc, -DBNAME => 'human');
  my $p_m_dbc = Bio::EnsEMBL::DBSQL::ProxyDBConnection->new(-DBC => $dbc, -DBNAME => 'mouse');
  
  # With a 10 minute timeout reconnection in milliseconds
  my $p_h_rc_dbc = Bio::EnsEMBL::DBSQL::ProxyDBConnection->new(-DBC => $dbc, -DBNAME => 'human', -RECONNECT_INTERVAL => (10*60*1000));

=head1 DESCRIPTION

This class is used to maintain one active connection to a database whilst it
appears to be working against multiple schemas. It does this by checking the
currently connected database before performing any query which could require
a database change such as prepare.

This class is only intended for internal use so please do not use unless
you are aware of what it will do and what the consequences of its usage are.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::ProxyDBConnection;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Utils::Proxy/;

use Bio::EnsEMBL::Utils::Argument qw/rearrange/;
use Bio::EnsEMBL::Utils::Exception qw/warning throw/;
use Bio::EnsEMBL::Utils::SqlHelper;

use Time::HiRes qw/time/;

sub new {
  my ($class, @args) = @_;
  my ($dbc, $dbname, $reconnect_interval) = rearrange([qw/DBC DBNAME RECONNECT_INTERVAL/], @args);
  throw "No DBConnection -DBC given" unless $dbc;
  throw "No database name -DBNAME given" unless $dbname;
  my $self = $class->SUPER::new($dbc);
  $self->dbname($dbname);
  if($reconnect_interval) {
    $self->reconnect_interval($reconnect_interval);
    $self->_last_used();
  }
  return $self;  
}

=head2 switch_database

  Description : Performs a switch of the backing DBConnection if the currently
                connected database is not the same as the database this proxy
                wants to connect to. It currently supports MySQL, Oracle and
                Postges switches is untested with all bar MySQL. If it
                cannot do a live DB/schema switch then it will disconnect
                the connection and then wait for the next process to
                connect therefore switching the DB.
  Exceptions  : None but will warn if you attempt to switch a DB with 
                active kids attached to the proxied database handle.

=cut

sub switch_database {
  my ($self) = @_;
  my $proxy = $self->__proxy();
  my $backing_dbname = $proxy->dbname();
  my $dbname = $self->dbname();
  
  my $switch = 0;
  if(defined $dbname) {
    if(defined $backing_dbname) {
      $switch = ($dbname ne $backing_dbname) ? 1 : 0;
    }
    else {
      $switch = 1;
    }
  }
  else {
    $switch = 1 if defined $backing_dbname;
  }
  
  if($switch) {
    $proxy->dbname($dbname);
    if($proxy->connected()) {
      my $kids = $proxy->db_handle()->{Kids};
      my $driver = lc($proxy->driver());
      #Edit to add other DB switching strategies on a per driver basis
      if($driver eq 'mysql') {
        $proxy->do('use '.$dbname);
      }
      elsif($driver eq 'oracle') {
        $proxy->do('ALTER SESSION SET CURRENT_SCHEMA = '.$dbname);
      }
      elsif($driver eq 'pg') {
        $proxy->do('set search_path to '.$dbname);
      }
      else {
        if($kids > 0) {
          warning "Attempting a database switch from '$backing_dbname' to '$dbname' with $kids active handle(s). Check your logic or do not use a ProxyDBConnection";
        }
        $proxy->disconnect_if_idle();
      }
    }
  }
  
  return $switch;
}

=head2 check_reconnection

  Description : Looks to see if the last time we used the backing DBI
                connection was greater than the reconnect_interval()
                provided at construction or runtime. If enought time has 
                elapsed then a reconnection is attempted. We do not
                attempt a reconnection if:
                
                  - No reconnect_interval was set
                  - The connection was not active
                
  Exceptions  : None apart from those raised from the reconnect() method
                from DBConnection
=cut

sub check_reconnection {
  my ($self) = @_;
  #Return early if we had no reconnection interval
  return unless $self->{reconnect_interval};
  
  my $proxy = $self->__proxy();
  
  #Only attempt it if we were connected; otherwise we can just skip
  if($proxy->connected()) {
    if($self->_require_reconnect()) {   
      $proxy->reconnect();
    }
    $self->_last_used();
  }
  return;
}

# Each time this is called we record the current time in seconds
# to be used by the _require_reconnect() method
sub _last_used {
  my ($self) = @_;
  $self->{_last_used} = int(time()*1000);
  return;
}

# Uses the _last_used() time and the current reconnect_interval() to decide
# if the connection has been unused for long enough that we should attempt
# a reconnect
sub _require_reconnect {
  my ($self) = @_;
  my $interval = $self->reconnect_interval();
  return unless $interval;
  my $last_used = $self->{_last_used};
  my $time_elapsed = int(time()*1000) - $last_used;
  return $time_elapsed > $interval ? 1 : 0;
}

=head2 reconnect_interval

	Arg[1]      : Integer reconnection interval in milliseconds
  Description : Accessor for the reconnection interval expressed in milliseconds
  Returntype  : Int miliseconds for a reconnection interval

=cut

sub reconnect_interval {
  my ($self, $reconnect_interval) = @_;
  $self->{'reconnect_interval'} = $reconnect_interval if defined $reconnect_interval;
  return $self->{'reconnect_interval'};
}

=head2 dbname

	Arg[1]      : String DB name
  Description : Accessor for the name of the database we should use whenever
                a DBConnection request is made via this class
  Returntype  : String the name of the database which we should use
  Exceptions  : None

=cut

sub dbname {
  my ($self, $dbname) = @_;
  $self->{'dbname'} = $dbname if defined $dbname;
  return $self->{'dbname'};
}

my %SWITCH_METHODS = map { $_ => 1 } qw/
  connect
  db_handle
  do
  prepare
  reconnect
  work_with_db_handle
/;


# Manual override of the SqlHelper accessor to ensure it always gets the Proxy
sub sql_helper {
  my ($self) = @_;
  if(! exists $self->{_sql_helper}) {
    my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(-DB_CONNECTION => $self);
    $self->{_sql_helper} = $helper;
  }
  return $self->{_sql_helper};
}

sub __resolver {
  my ($self, $package, $method) = @_;
  if($self->__proxy()->can($method)) {
    if($SWITCH_METHODS{$method}) {
      return sub {
        my ($local_self, @args) = @_;
        $local_self->check_reconnection();
        $local_self->switch_database();
        $local_self->_last_used();
        return $local_self->__proxy()->$method(@args);
      };
    }
    else {
      return sub {
        my ($local_self, @args) = @_;
        return $local_self->__proxy()->$method(@args);
      };
    }
  }
  return;
}

1;

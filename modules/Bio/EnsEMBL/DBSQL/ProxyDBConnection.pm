=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::ProxyDBConnection - Database connection wrapper allowing
for one backing connection to be used for multiple DBs

=head1 SYNOPSIS

  my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(-HOST => 'host', -PORT => 3306, -USER => 'user');
  my $p_h_dbc = Bio::EnsEMBL::DBSQL::ProxyDBConnection->new(-DBC => $dbc, -DBNAME => 'human');
  my $p_m_dbc = Bio::EnsEMBL::DBSQL::ProxyDBConnection->new(-DBC => $dbc, -DBNAME => 'mouse');

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

use Bio::EnsEMBL::Utils::Exception qw/warning/;
use Bio::EnsEMBL::Utils::Argument qw/rearrange/;

use base qw/Bio::EnsEMBL::Utils::Proxy/;

sub new {
  my ($class, @args) = @_;
  my ($dbc, $dbname) = rearrange([qw/DBC DBNAME/], @args);
  my $self = $class->SUPER::new($dbc);
  $self->dbname($dbname);
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

sub __resolver {
  my ($self, $package, $method) = @_;
  if($self->__proxy()->can($method)) {
    if($SWITCH_METHODS{$method}) {
      return sub {
        my ($local_self, @args) = @_;
        $local_self->switch_database();
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

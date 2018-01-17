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

package XrefMapper::ChecksumMapper;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::DBSQL::DBConnection;

use base qw(XrefMapper::BasicMapper);

my $DEFAULT_METHOD = 'XrefMapper::Methods::MySQLChecksum';

sub new {
  my($class, $mapper) = @_;
  my $self = bless {}, $class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->mapper($mapper);
  return $self;
}

sub _xref_helper {
  my ($self) = @_;
  return $self->xref()->dbc()->sql_helper();
}

sub logic_name {
  my ($self) = @_;
  return 'xrefchecksum';
}

sub mapper {
  my ($self, $mapper) = @_;
  $self->{mapper} = $mapper if defined $mapper;
  return $self->{mapper};
}

sub method {
  my ($self, $method) = @_;
  $self->{method} = $method if defined $method;
  return $self->{method};
}

sub verbose {
  my ($self) = @_;
  return $self->mapper()->verbose();
}

# No default target file, implemented by subclasses where necessary
sub target {
  return;
}

sub process {
  my ($self, $db_url, $species_id) = @_;

  $self->_update_status('checksum_xrefs_started');
  my $source_id = $self->source_id();
  my $target = $self->target();
  my $object_type = $self->object_type;

  if($self->_map_checksums($db_url)) {
    my $method = $self->get_method();
    my $results = $method->run($target, $source_id, $object_type, $db_url);
    $self->log_progress('Starting upload');
    $self->upload($results, $species_id);
  }

  $self->_update_status('checksum_xrefs_finished');
  return;
}

sub upload {
  my ($self, $results, $species_id) = @_;
  #The elements come in as an array looking like
  #  [ { id => 1, upi => 'UPI00000A', object_type => 'Translation' } ]
  
  my $insert_xref = <<'SQL';
INSERT INTO xref (source_id, accession, label, version, species_id, info_type)
values (?,?,?,?,?,?)
SQL
  my $insert_object_xref = <<'SQL';
INSERT INTO object_xref (ensembl_id, ensembl_object_type, xref_id, linkage_type, ox_status)
values (?,?,?,?,?)
SQL
  
  my $h = $self->_xref_helper();
  my $source_id = $self->source_id();
  $species_id = $self->species_id() unless defined $species_id;
  if (!defined $species_id) { return; }
   
  $h->transaction(-CALLBACK => sub {
    
    $self->log_progress('Deleting records from previous possible upload runs');
    $self->_delete_entries('object_xref');
    $self->_delete_entries('xref');
    
    $self->log_progress('Starting xref insertion');
    #Record UPIs to make sure we do not attempt to insert duplicate UPIs
    my %upi_xref_id;
    $h->batch(-SQL => $insert_xref, -CALLBACK => sub {
      my ($sth) = @_;
      foreach my $e (@{$results}) {
        my $upi = $e->{upi};
        if(exists $upi_xref_id{$upi}) {
          $e->{xref_id} = $upi_xref_id{$upi};
        }
        else {
          $sth->execute($source_id, $e->{upi}, $e->{upi}, 1, $species_id, 'CHECKSUM');
          my $id = $sth->{'mysql_insertid'};
          $e->{xref_id} = $id;
          $upi_xref_id{$upi} = $id;
        }
      }
      return;
    });
    
    $self->log_progress('Starting object_xref insertion');
    $h->batch(-SQL => $insert_object_xref, -CALLBACK => sub {
      my ($sth) = @_;
      foreach my $e (@{$results}) {
        $sth->execute($e->{id}, $e->{object_type}, $e->{xref_id}, 'CHECKSUM', 'DUMP_OUT');
      }
      return;
    });
  });
  
  $self->log_progress('Finished insertions');
  
  return;
}

sub _delete_entries {
  my ($self, $table) = @_;
  $self->log_progress('Deleting entries from %s', $table);
  my $lookup = {
    xref => <<'SQL',
DELETE  x
FROM    xref x
WHERE   x.source_id    = ?
SQL
    object_xref => <<'SQL',
DELETE  ox
FROM    xref x,
        object_xref ox
WHERE   x.source_id    = ?
AND     ox.xref_id          = x.xref_id
SQL
  };
  
  my $sql = $lookup->{$table};
  throw "Cannot find delete SQL for the table $table" unless $sql;
  my $source_id = $self->source_id();
  my $count = $self->_xref_helper()->execute_update(-SQL => $sql, -PARAMS => [$source_id]);
  my $type = ($count == 1) ? 'entry' : 'entries';
  $self->log_progress('Deleted %s %s from %s', $count, $type, $table);
  return;
}

sub source_id {
  my ($self) = @_;
  return $self->_xref_helper()->execute_single_result(
    -SQL => 'select source_id from source where name=?',
    -PARAMS => [$self->external_db_name()]
  );
}

sub species_id {
  my ($self) = @_;
  my $species_id = $self->SUPER::species_id();
  if(! defined $species_id) {
    $species_id = $self->get_id_from_species_name($self->core()->species());
    $self->SUPER::species_id($species_id);
  }
  return $species_id;
}

sub get_method {
  my ($self) = @_;
  my $method_class = $DEFAULT_METHOD;
  eval "require ${method_class};";
  if($@) {
    throw "Cannot require the class ${method_class}. Make sure your PERL5LIB is correct: $@";
  }
  return $method_class->new( -MAPPER => $self );
}

############# INTERNAL METHODS

sub _update_status {
  my ($self, $status) = @_;
  if($self->xref()) {
    my $h = $self->_xref_helper();
    my $sql = q{insert into process_status (status, date) values(?,now())};
    $h->execute_update(-SQL => $sql, -PARAMS => [$status]);
  }
  else {
    my $time = localtime();
    $self->log_progress(q{Status Update '%s' @ %s}."\n", $status, $time);
  }
  return;
}

sub _map_checksums {
  my ($self, $db_url) = @_;
  my $source_id = $self->source_id();
  my $dbc = $self->mapper->xref->dbc;
  if (defined $db_url) { 
    $source_id = 1;
    my ($dbconn_part, $driver, $user, $pass, $host, $port, $dbname, $table_name, $tparam_name, $tparam_value, $conn_param_string) =
            $db_url =~ m{^((\w*)://(?:(\w+)(?:\:([^/\@]*))?\@)?(?:([\w\-\.]+)(?:\:(\d*))?)?/([\w\-\.]*))(?:/(\w+)(?:\?(\w+)=([\w\[\]\{\}]*))?)?((?:;(\w+)=(\w+))*)$};
    $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
      -dbname => $dbname,
      -user => $user,
      -pass => $pass,
      -host => $host,
      -port => $port);
  }
  my $count = $dbc->sql_helper()->execute_single_result(-SQL => 'select count(*) from checksum_xref where source_id = ' . $source_id);
  return $count;
}

sub log_progress {
  my ( $self, $fmt, @params ) = @_;
  return if (!$self->verbose);
  printf( STDERR "CHKSM==> %s\n", sprintf( $fmt, @params ) );
}

1;

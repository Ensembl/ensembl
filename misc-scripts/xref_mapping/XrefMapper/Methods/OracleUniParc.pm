package XrefMapper::Methods::OracleUniParc;

use strict;
use warnings;

use base qw/XrefMapper::Methods::ChecksumBasic/;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::SqlHelper;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use List::Util qw(max);

my $DEFAULT_BATCH_SIZE = 10000;

my $UNIPARC_SQL = <<'SQL';
SELECT p.UPI
FROM UNIPARC.PROTEIN p 
WHERE p.md5 = ?
SQL

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($batch_size) = rearrange([qw(batch_size)], @args);
  if(! $batch_size) {
    $self->batch_size($DEFAULT_BATCH_SIZE);
  }
  return $self;
}

sub checksum {
  my ($self, $sequence) = @_;
  return uc($self->md5_checksum($sequence));
}

sub perform_mapping {
  my ($self, $sequences) = @_;
  
  my @final_results;
  
  $self->oracle_dbc()->sql_helper()->batch(-SQL => $UNIPARC_SQL, -CALLBACK => sub {
    my ($sth) = @_;
    foreach my $sequence (@{$sequences}) {
      my $checksum = $self->checksum($sequence);
      $sth->execute($checksum);
      my $upi;
      while(my $row = $sth->fetchrow_arrayref()) {
        my ($local_upi) = @{$row};
        if(defined $upi) {
          throw sprintf('The sequence %s had a checksum of %s but this resulted in more than one UPI: [%s, %s]', $sequence->id(), $checksum, $upi, $local_upi);
        }
        $upi = $local_upi;
      }
      if(defined $upi){
        push(@final_results, { id => $sequence->id(), upi => $upi, object_type => 'Translation' });
      }
    }
    return;
  });
  
  return \@final_results;
}

sub oracle_dbc {
  my ($self) = @_;
  if(! exists $self->{oracle_dbc}) {
    my $dbc = $self->mapper()->uniparc()->dbc();
    $dbc->disconnect_when_inactive(0);
    $dbc->driver('Oracle');
    $self->{oracle_dbc} = $dbc;
  }
  return $self->{oracle_dbc};
}

sub DESTROY {
  my ($self) = @_;
  $self->oracle_dbc()->disconnect_if_idle() if $self->oracle_dbc();
  return;
}

1;
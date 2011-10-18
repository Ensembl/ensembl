package Bio::EnsEMBL::DBSQL::DataFileAdaptor;

=pod

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

Bio::EnsEMBL::DBSQL::DataFileAdaptor

=head1 SYNOPSIS

	my $dfa = $dba->get_DataFileAdaptor();
	my $file = $dfa->fetch_by_dbID(1);
	my $files = $dfa->fetch_all();

=head1 DESCRIPTION

Provides a database wrapper to store the locations of files and to pull these
records back out. DataFile objects can only provide basic information but they
can return an intended external database adaptor which can be used to 
parse the information. This system assumes nothing about the file just that
your parser can access it.

Files are supported over any protocol your parser supports and locations can be
made absolute, built on the fly or versioned.

=head1 METHODS

=cut

use strict;
use warnings;

use base qw/Bio::EnsEMBL::DBSQL::BaseAdaptor/;

use Bio::EnsEMBL::DataFile;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw/throw warning/;
use Bio::EnsEMBL::Utils::Scalar qw/:assert/;

=head2 DataFile_to_extension

  Arg[1]      : Bio::EnsEMBL::DataFile
  Example     : my $ext = $dfa->DataFile_to_extension($bam_df);
  Description : Returns an expected extension for the given DataFile type
  Returntype  : Scalar of the expected file extension
  Exceptions  : Raised if the given file type is not understood

=cut

sub DataFile_to_extension {
  my ($self, $df) = @_;
  my $type = $df->file_type();
  my $ext = {
    BAM     => 'bam',
#    BIGBED  => 'bb',
    BIGWIG  => 'bw',
    VCF     => 'vcf',
  }->{$type}; 
  throw sprintf(q{No extension found for the type '%s'}, $type ) if ! $ext;
  return $ext;
}

=head2 DataFile_to_adaptor

  Arg[1]     	: Bio::EnsEMBL::DataFile
  Arg[2]      : (optional) base path
  Example			: my $bam = $dfa->DataFile_to_adaptor($bam_df);
  Description	: Returns an adaptor instance which will access the given DataFile
  Returntype 	: Scalar actual return depends upon the given file type
  Exceptions 	: Raised if the given file type is not understood

=cut

sub DataFile_to_adaptor {
  my ($self, $df, $base) = @_;
  my $type = $df->file_type();
  my $dispatch = {
    BAM     => sub {
      require Bio::EnsEMBL::ExternalData::BAM::BamAdaptor;
      return Bio::EnsEMBL::ExternalData::BAM::BamAdaptor->new($df->path($base));
    },
#    BIGBED  => sub {
#      require Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor;
#      return Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor->new($df->path($base));
#    },
    BIGWIG  => sub {
      require Bio::EnsEMBL::ExternalData::BigFile::BigWigAdaptor;
      return Bio::EnsEMBL::ExternalData::BigFile::BigWigAdaptor->new($df->path($base));
    },
    VCF     => sub {
      require Bio::EnsEMBL::ExternalData::VCF::VCFAdaptor;
      return Bio::EnsEMBL::ExternalData::VCF::VCFAdaptor->new($df->path($base));
    },
  }->{$type}; 
  throw sprintf(q{No handler found for the type '%s'}, $type ) if ! $dispatch;
  return $dispatch->();
}


sub fetch_all_by_Analysis {
  my ($self, $analysis) = @_;
  assert_ref($analysis, 'Bio::EnsEMBL::Analysis', 'analysis');
  $self->bind_param_generic_fetch($analysis->dbID(), SQL_INTEGER);
  return $self->generic_fetch('df.analysis_id =?');
}

sub fetch_all_by_CoordSystem {
  my ($self, $cs) = @_;
  assert_ref($cs, 'Bio::EnsEMBL::CoordSystem', 'coord_system');
  $self->bind_param_generic_fetch($cs->dbID(), SQL_INTEGER);
  return $self->generic_fetch('df.coord_system_id =?');
}


sub fetch_by_name_and_type {
  my ($self, $name, $type) = @_;
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  $self->bind_param_generic_fetch($type, SQL_VARCHAR);
  my $results = $self->generic_fetch('df.name =? and df.file_type =?');
  return $results->[0] if @{$results};
  return;
}

sub generic_fetch {
  my ($self, $constraint) = @_;
  $constraint ||= q{};
  
  my $sql = <<'SQL';
select df.data_file_id, df.coord_system_id, df.analysis_id, df.name, df.version_lock, df.absolute, df.url, df.file_type
from data_file df
join coord_system cs using (coord_system_id) 
where cs.species_id =?
SQL
  $sql .= 'AND '.$constraint if $constraint;
  
  my $params = $self->bind_param_generic_fetch();
  if(defined $params) {
    $self->{'_bind_param_generic_fetch'} = ();
  }
  else {
    $params = [];
  }
  unshift(@{$params}, $self->db()->species_id());
  
  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $aa = $self->db()->get_AnalysisAdaptor();
  
  return $self->dbc()->sql_helper()->execute(-SQL => $sql, -PARAMS => $params, -CALLBACK => sub {
    my ($row) = @_;
    my ($data_file_id, $coord_system_id, $analysis_id, $name, $version_lock, $absolute, $url, $file_type) = @{$row};
    my $hash = {
      dbID          => $data_file_id,
      adaptor       => $self,
      coord_system  => $csa->fetch_by_dbID($coord_system_id),
      analysis      => $aa->fetch_by_dbID($analysis_id),
      name          => $name,
      version_lock  => $version_lock,
      absolute      => $absolute,
      file_type     => $file_type,
    };
    $hash->{url} = $url if $url;
    return Bio::EnsEMBL::DataFile->new_fast($hash);
  });
}

sub store {
  my ($self, $df) = @_;
  
  assert_ref($df, 'Bio::EnsEMBL::DataFile', 'datafile');
  
  if ($df->is_stored($self->db())) {
    return $df->dbID();
  }
  
  throw 'Analysis is not defined for this data file' if ! defined $df->analysis();
  throw 'Coord system is not defined for this data file' if ! defined $df->coord_system();
  
  my $sql = <<'SQL';
INSERT INTO data_file (coord_system_id, analysis_id, name, version_lock, absolute, url, file_type)
VALUES (?,?,?,?,?,?,?)
SQL
  my $params = [
    [$df->coord_system()->dbID(), SQL_INTEGER], 
    [$df->analysis()->dbID(), SQL_INTEGER],
    [$df->name(), SQL_VARCHAR],
    [$df->version_lock(), SQL_INTEGER],
    [$df->absolute(), SQL_INTEGER],
    [$df->url(), SQL_VARCHAR],
    [$df->file_type(), SQL_VARCHAR],
  ];
  $self->dbc()->sql_helper()->execute_update(-SQL => $sql, -PARAMS => $params, -CALLBACK => sub {
    my ( $sth, $dbh ) = @_;
    $df->dbID($self->last_insert_id());
    return;
  });
  $df->adaptor($self);
  
  return $df->dbID();
}

sub update {
  my ($self, $df) = @_;
  
  assert_ref($df, 'Bio::EnsEMBL::DataFile', 'datafile');
  
  if (! $df->is_stored($self->db())) {
    $self->store($df);
    return;
  }
  
  my $sql = <<'SQL';
UPDATE data_file SET coord_system_id =?, analysis_id=?, name=?, version_lock=?, absolute=?, url=?, file_type=?
WHERE data_file_id =?
SQL
  my $params = [
    [$df->coord_system()->dbID(), SQL_INTEGER], 
    [$df->analysis()->dbID(), SQL_INTEGER],
    [$df->name(), SQL_VARCHAR],
    [$df->version_lock(), SQL_INTEGER],
    [$df->absolute(), SQL_INTEGER],
    [$df->url(), SQL_VARCHAR],
    [$df->file_type(), SQL_VARCHAR],
    [$df->dbID(), SQL_INTEGER],
  ];
  $self->dbc()->sql_helper()->execute_update(-SQL => $sql, -PARAMS => $params);
  return;
}

sub delete {
  my ($self, $df) = @_;
  
  assert_ref($df, 'Bio::EnsEMBL::DataFile', 'datafile');
  
  if (! $df->is_stored($self->db())) {
    throw "Cannot delete the data file if it has not already been stored in this database";
  }
  
  $self->dbc()->sql_helper()->execute_update(
    -SQL => 'DELETE from data_file where data_file_id =?', 
    -PARAMS => [[$df->dbID(), SQL_INTEGER]],
  );
  
  return;
}

sub _tables {
  my ($self) = @_;
  return (
    [qw/data_file df/]
  );
}

1;
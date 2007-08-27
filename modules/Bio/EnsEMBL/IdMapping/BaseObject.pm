package Bio::EnsEMBL::IdMapping::BaseObject;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);


=head2 new

  Arg[1]      : 
  Example     : 
  Description : constructor
  Return type : 
  Exceptions  : 
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($logger, $conf, $cache) = rearrange(['LOGGER', 'CONF', 'CACHE'], @_);

  unless ($logger->isa('Bio::EnsEMBL::Utils::Logger')) {
    throw("You must provide a Bio::EnsEMBL::Utils::Logger for logging.");
  }
  
  unless ($conf->isa('Bio::EnsEMBL::Utils::ConfParser')) {
    throw("You must provide configuration as a Bio::EnsEMBL::Utils::ConfParser object.");
  }
  
  unless ($cache->isa('Bio::EnsEMBL::IdMapping::Cache')) {
    throw("You must provide configuration as a Bio::EnsEMBL::IdMapping::Cache object.");
  }
  
  my $self = {};
  bless ($self, $class);

  # initialise
  $self->logger($logger);
  $self->conf($conf);
  $self->cache($cache);
  
  return $self;
}


sub get_filehandle {
  my $self = shift;
  my $filename = shift;
  my $path_append = shift;
  my $mode = shift;

  throw("Need a filename for this filehandle.") unless (defined($filename));
  
  my $path = $self->conf->param('dumppath');
  $path = path_append($path, $path_append) if (defined($path_append));

  $mode ||= '>';
  
  open(my $fh, $mode, "$path/$filename") or
    throw("Unable to open $path/$filename: $!");

  return $fh;
}


sub file_exists {
  my $self = shift;
  my $filename = shift;
  my $path_append = shift;

  my $path = $self->conf->param('dumppath');
  $path = path_append($path, $path_append) if (defined($path_append));

  return (-s "$path/$filename");
}


sub fetch_value_from_db {
  my $self = shift;
  my $dbh = shift;
  my $sql = shift;

  throw("Need a db handle.") unless ($dbh and $dbh->isa('DBI::db'));
  throw("Need an SQL query to execute.") unless ($sql);

  my $sth = $dbh->prepare($sql);
  $sth->execute;
  my ($retval) = $sth->fetchrow_array;

  return $retval;
}


sub dump_table_to_file {
  my $self = shift;
  my $dbtype = shift;
  my $table = shift;
  my $filename = shift;

  # argument check
  unless (($dbtype eq 'source') or ($dbtype eq 'target')) {
    throw("Missing or unknown db type: $dbtype.");
  }
  throw("Need a table name.") unless ($table);
  throw("Need a filename.") unless ($filename);
  
  my $fh = $self->get_filehandle($filename, 'tables');

  my $dba = $self->cache->get_DBAdaptor($dbtype);
  my $dbh = $dba->dbc->db_handle;
  my $sth = $dbh->prepare("SELECT * FROM $table");
  $sth->execute;

  my $i = 0;

  while (my @row = $sth->fetchrow_array) {
    $i++;

    # use '\N' for NULL values
    for (my $j = 0; $j < scalar(@row); $j++) {
      $row[$j] = '\N' unless (defined($row[$j]));
    }
    
    print $fh join("\t", @row);
    print $fh "\n";
  }

  $sth->finish;
  
  return $i;
}


sub upload_file_into_table {
  my $self = shift;
  my $dbtype = shift;
  my $table = shift;
  my $filename = shift;

  # argument check
  unless (($dbtype eq 'source') or ($dbtype eq 'target')) {
    throw("Missing or unknown db type: $dbtype.");
  }
  throw("Need a table name.") unless ($table);
  throw("Need a filename.") unless ($filename);

  # sanity check for dry run
  if ($self->conf->param('dry_run')) {
    $self->logger->warning("dry_run - skipping db upload for $filename.\n");
    return;
  }
  
  my $file = join('/', $self->conf->param('dumppath'), 'tables', $filename);
  
  if (-s $file) {
    
    my $dba = $self->cache->get_DBAdaptor($dbtype);
    my $dbh = $dba->dbc->db_handle;
    my $sql = qq(LOAD DATA LOCAL INFILE '$file' INTO TABLE $table);
    my $sth = $dbh->prepare($sql);
    $sth->execute;
    $sth->finish;

  } else {
    $self->logger->warning("No data found in file $filename.\n", 1);
  }

}


sub logger {
  my $self = shift;
  $self->{'_logger'} = shift if (@_);
  return $self->{'_logger'};
}


sub conf {
  my $self = shift;
  $self->{'_conf'} = shift if (@_);
  return $self->{'_conf'};
}


sub cache {
  my $self = shift;
  $self->{'_cache'} = shift if (@_);
  return $self->{'_cache'};
}


1;


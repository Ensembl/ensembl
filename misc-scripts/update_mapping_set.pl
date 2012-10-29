#!/usr/local/ensembl/bin/perl

=head1 NAME

update_mapping_set.pl   - script to update 2 tables: mapping_set and
                          seq_region_id mapping.  These tables are
                          supposed to store different mapping of
                          seq_region_id between releases

=head1 SYNOPSIS

update_mapping_set.pl [arguments]

Required arguments:

  --user=user             username for the database

  --pass=pass             password for the database

  --release=num_release   release we want to run the updating

Optional arguments:

  --host=host             server where the core databases are stored
                          (default: ens-staging)

  --oldhost = oldhost   server where the old release databases are stored

  --dbname=dbname         if you want a single database to update
                          the mapping_set information (default: all
                          databases)

  --port=port             port (default: 3306)

  --oldport=port          old database server port (default: 3306)

  --olduser=user          old database server username

  --oldpass=pass          password for old database server

  --help                  print help (this message)

=head1 DESCRIPTION

This script will update the mapping_set table with the current
schema_build information and indicate if the seq_region table has
changed from previous release.

If it has, the table seq_region_mapping will contain the
relation between old (external_seq_region_id) and the current
(internal_seq_region_id).

=head1 EXAMPLES

Update mapping_set information for all databases in ens-staging in
release NN (the usual use case in release process):

  $ ./update_mapping_set.pl --user ensadmin \
    --pass password --release NN --old_host ensembldb-ensembl.org

Update mapping_set information only for pig database in ens-genomics1:

  $ ./update_mapping_set.pl --host ens-genomics1 \
    --user ensadmin --pass password --dbname my_pig_db --release 52

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Daniel Rios <dani@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use DBI qw(:sql_types);

use Bio::EnsEMBL::Utils::Exception qw(throw);

## Command line options

my $host    = 'ens-staging';
my $host2   = 'ens-staging2';
my $oldhost = 'ens-livemirror';
my $dbname = undef;
my $user = undef;
my $pass = undef;
my $port = 3306;
my $oldport = 3306;
my $olduser = "ensro";
my $oldpass = undef;
my $previous_dbname = undef;
my $help = undef;
my $release = undef;
my $dry_run = undef;

GetOptions('host=s'    => \$host,
	   'dbname=s'  => \$dbname,
	   'user=s'    => \$user,
	   'pass=s'    => \$pass,
	   'port=s'    => \$port,
           'release=i' => \$release,
	   'oldhost=s' => \$oldhost,
	   'oldport=s' => \$oldport,
	   'olduser=s' => \$olduser,
	   'oldpass=s' => \$oldpass,
           'previous_dbname=s'  => \$previous_dbname,
	   'help'      => \$help,
	   'dry_run'   => \$dry_run,
	   );

pod2usage(1) if($help);
throw("--user argument required") if (!defined($user));
throw("--release argument required") if(!defined($release));

my $database = 'information_schema';
my $old_dbh =  DBI->connect("DBI:mysql:database=$database;host=$oldhost;port=$oldport",$olduser,$oldpass);

foreach my $h ($host,$host2) {
  my $dbh = DBI->connect("DBI:mysql:database=$database;host=$h;port=$port",$user,$pass);

  #since there is no database defined, will run it agains all core databases
  my $pattern;
  if (!defined ($dbname)){
    $pattern = "_core_".$release."_";
  }
  else{
    $pattern = $dbname;
  }

# Fetch all databases matching the pattern
  my $sth = $dbh->prepare("SHOW DATABASES like \'%$pattern%\'");
  $sth->execute();
  my $dbs = $sth->fetchall_arrayref();
  foreach my $db_name (@{$dbs}){
    my $current_dbname = $db_name->[0];
    print STDERR "Going to update mapping for $current_dbname....\n";
    my $mapping_set_id = get_max_mapping_set_id($dbh,$current_dbname) + 1;
    my $sth_seq_mapping = $dbh->prepare("INSERT INTO $current_dbname.seq_region_mapping VALUES(?,?,?)");
    my $sth_mapping_set = $dbh->prepare("INSERT INTO $current_dbname.mapping_set VALUES(?,?,?)");
    my $sth_update_build = $dbh->prepare("UPDATE $current_dbname.mapping_set SET internal_schema_build = ?");
    my $sth_update_old = $dbh->prepare("UPDATE $current_dbname.seq_region_mapping SET internal_seq_region_id = ? WHERE internal_seq_region_id = ?");
    my $sth_remove_deprecated = $dbh->prepare("DELETE FROM $current_dbname.seq_region_mapping WHERE internal_seq_region_id = ?");
    my $schema_build = get_schema_and_build($current_dbname);
    my $current_assembly = get_assembly($dbh,$current_dbname) ;
    my $count_removed = 0;
    my $count_updated = 0;
    my $count_added = 0;

    $sth_update_build->execute($schema_build) unless $dry_run;
    if (!$previous_dbname) {
       $previous_dbname = &get_previous_dbname($old_dbh,$current_dbname,$release);
    }

# If there is no previous database, no mapping needed
    if (!defined($previous_dbname)) {
       print STDERR "First instance known for $current_dbname, no mapping needed\n";
       next;
    }

# If it is a new assembly, no mapping needed
    my $old_assembly = get_assembly($old_dbh,$previous_dbname); 
    if ($old_assembly ne $current_assembly) { 
      print STDERR "New assembly $current_assembly for $current_dbname, no mapping needed\n" ;
      next; 
    }

    my $previous_schema_build = get_schema_and_build($previous_dbname);
    my $new_mapping = $sth_mapping_set->execute($mapping_set_id,$schema_build,$previous_schema_build) unless $dry_run;

    if (!$new_mapping) {
      print STDERR "Mapping already run for this schema_build, please remove entry before proceeding\n" ;
      exit;
    }

# If there has been no change in seq_region, no mapping needed
    my $cur_seq_region_checksum = &get_seq_region_checksum($dbh,$current_dbname);
    my $previous_seq_region_checksum = &get_seq_region_checksum($old_dbh,$previous_dbname);
    if ($cur_seq_region_checksum == $previous_seq_region_checksum) { 
      print STDERR "No change in seq_region for $current_dbname, no mapping needed\n";
      next; 
    }

# There has been a seq_region change between releases, add the relation old_seq_region_id->new_seq_region_id
    my $current_seq_region =  &read_seq_region($dbh,$current_dbname);
    my $old_seq_region = &read_seq_region($old_dbh,$previous_dbname);

# Update the seq_region_mapping table with the old->new seq_region_id relation
    foreach my $seq_region_name (keys %{$old_seq_region}){
      my $current_name_hash = $current_seq_region->{$seq_region_name};
      my $old_name_hash = $old_seq_region->{$seq_region_name};

# If the seq_region has disappeared, remove previous entries for that id
      if (!defined $current_name_hash) {
        my $id = get_seq_region_id($old_dbh,$previous_dbname, $seq_region_name);
        $count_removed += $sth_remove_deprecated->execute($id) unless $dry_run;
        next;
      }
      foreach my $length (keys %{$old_name_hash}){
        my $current_length_hash = $current_name_hash->{$length};
        my $old_length_hash = $old_name_hash->{$length};

# The seq_region might have a different length
        if (!defined $current_length_hash) {
          next;
        }
        foreach my $cs (keys %{$old_length_hash}) {
          my $current_cs_hash = $current_length_hash->{$cs};
          my $old_cs_hash = $old_length_hash->{$cs};

# The coord system might have changed
          if (!defined $current_cs_hash) {
            next;
          }
          foreach my $id (keys %{$old_cs_hash}) {
            my $current_id = $current_cs_hash->{$id};
            my $old_id = $old_cs_hash->{$id};

# If no change, no need to write out
            if (!defined $current_id || $old_id == $current_id) {
              next;
            }

# If there is a change, update any existing entries for this seq_region to the new id
# Then, add a new entry to map said id to the old release
            $count_updated += $sth_update_old->execute($current_id,$old_id) unless $dry_run;
            $count_added += $sth_seq_mapping->execute($old_id,$current_id, $mapping_set_id) unless $dry_run;
          }
        }
      }
    }
    print STDERR "For $current_dbname, removed $count_removed, added $count_added, updated $count_updated seq_region_mapping entries\n\n" ;
  }
}


# For a given database, will return the seq_region_name->seq_region_id relation
sub read_seq_region {
  my ($dbh, $dbname) = @_;
  my (%seq_region_hash, $seq_region_id, $seq_region_name, $coord_system_id, $length, $cs_name, $cs_rank);
  my $sth = $dbh->prepare("SELECT seq_region_id, s.name, length, cs.name, cs.rank FROM $dbname.seq_region s, $dbname.coord_system cs WHERE cs.coord_system_id = s.coord_system_id");
  $sth->execute();
  $sth->bind_col(1,\$seq_region_id);
  $sth->bind_col(2,\$seq_region_name);
  $sth->bind_col(3,\$length);
  $sth->bind_col(4,\$cs_name);
  $sth->bind_col(5,\$cs_rank);
  while ($sth->fetch){
    $seq_region_hash{$seq_region_name}{$length}{$cs_name}{$cs_rank} = $seq_region_id;
  }
  return \%seq_region_hash;
}

# For a given database, returns the size of the seq_region_table
sub get_seq_region_checksum {
    my ($dbh, $dbname) = @_;
    my $sth_status = $dbh->prepare("checksum table $dbname.seq_region") ;
    $sth_status->execute();
    my $table_status = $sth_status->fetchrow_array();
    return $table_status; #return the size of the table
}

sub get_seq_region_id {
    my ($dbh, $dbname, $seq_region_name) = @_;
    my $sth_region = $dbh->prepare("SELECT seq_region_id FROM $dbname.seq_region WHERE name = ?");
    $sth_region->execute($seq_region_name);
    my $seq_region_id = $sth_region->fetchrow_array();
    return $seq_region_id;
}

# Will return the max mapping_set_id being used in the mapping_set table
sub get_max_mapping_set_id {
    my ($dbh, $dbname) = @_;
    my $sth_mapping = $dbh->prepare("select max(mapping_set_id) from $dbname.mapping_set");
    $sth_mapping->execute();
    my ($max_mapping_set_id) = $sth_mapping->fetchrow_array();
    if (!defined $max_mapping_set_id) { return 0; }
    return $max_mapping_set_id;
}

# This method will return the name of the previous database to release for same species (assuming is present)
sub get_previous_dbname {
    my ($dbh, $dbname, $release) = @_;
    my $previous_dbname;
    $dbname =~ /(^([a-z]+_){2,3}core_)/;
    if (!$1) { throw("Database name $dbname is not in the right format"); }
    my $previous_release_name = $1 . (--$release);
    my $previous_sth = $dbh->prepare("show databases like \'%$previous_release_name%\'");
    $previous_sth->execute();
    ($previous_dbname) = $previous_sth->fetchrow_array() ;
    return $previous_dbname;
}

# For a standard ensembl database name, returns the release number and assembly
sub get_schema_and_build {
  my ($dbname) = @_;
  my @dbname = split/_/, $dbname;
  return join "_", $dbname[($#dbname -1)], $dbname[($#dbname)];
}

# Returns the assembly name for a given database
sub get_assembly {
  my ($dbh, $dbname) = @_;
  my $sth = $dbh->prepare("select meta_value from $dbname.meta where meta_key = 'assembly.default'");
  $sth->execute();
  return $sth->fetchrow();
}

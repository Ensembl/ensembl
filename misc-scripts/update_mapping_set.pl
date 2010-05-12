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

  --oldport=port          old database server port (default: 5306)

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

If it hasn't, a single entry in the mapping_set table with the same
mapping_set_id as the previous release will be stored.

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

use constant INITIAL_MAPPING => 1;
use constant SAME_MAPPING => 2;
use constant NEW_MAPPING => 3;

## Command line options

my $host = 'ens-staging';
my $oldhost = 'ensembldb.ensembl.org';
my $dbname = undef;
my $user = undef;
my $pass = undef;
my $port = 3306;
my $oldport = 5306;
my $olduser = "anonymous";
my $oldpass = undef;
my $help = undef;
my $release = undef;

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
	   'help'      => \$help
	   );

pod2usage(1) if($help);
throw("--user argument required") if (!defined($user));
throw("--pass argument required") if (!defined($pass));
throw("--release argument required") if(!defined($release));

my $database = 'information_schema';
my $dbh = DBI->connect("DBI:mysql:database=$database;host=$host;port=$port",$user,$pass);
my $old_dbh =  DBI->connect("DBI:mysql:database=$database;host=$oldhost;port=$oldport",$olduser,$oldpass);
my $status;
my $database_name;

#since there is no database defined, will run it agains all core databases
my $pattern;
if (!defined ($dbname)){
    $pattern = "_core_".$release."_";
}
else{
$pattern = $dbname;
}
#fetch all databases matching the pattern
print STDERR $pattern."\n";
my $sth = $dbh->prepare("SHOW DATABASES WHERE `database` REGEXP \'$pattern\'");
$sth->execute();
my $dbs = $sth->fetchall_arrayref();
my $schema_build;
foreach my $db_name (@{$dbs}){
  print STDERR "Going to update mapping for $db_name->[0]....\n";
  my $mapping_set_id;
  my $current_seq_region = (); # hash containing the relation seq_region_name->seq_region_id for the current database
  my $old_seq_region = (); #hash containing the previous database relation seq_region_name->seq_region_id
  my $sth_seq_mapping = $dbh->prepare("INSERT INTO $db_name->[0].seq_region_mapping VALUES(?,?,?)");
  my $sth_mapping_set = $dbh->prepare("INSERT INTO $db_name->[0].mapping_set VALUES(?,?)");
  $status = &mapping_status($dbh,$old_dbh, $db_name->[0],\$mapping_set_id,$release);
  $schema_build = get_schema_and_build($db_name->[0]);

  #add mapping_set information
  if ($status == INITIAL_MAPPING){
    $sth_mapping_set->execute($mapping_set_id,$schema_build);
    #first time run the script, create new entry in mapping_set and copy seq_region entries in seq_region_mapping

    ############
    #Actually NO only store the differences so for the initial one it is NONE.
    ############

    #        $current_seq_region =  &read_seq_region($dbh,$db_name->[0]);
    #        #copy the seq_region_id in the seq_region_mapping
    #        foreach my $seq_region_name (keys %{$current_seq_region}){
    #            #when copying there won't be any ambiguity with coord_systems
    #            foreach my $region_id (values %{$current_seq_region->{$seq_region_name}}){
    #                $sth_seq_mapping->execute($region_id,$region_id,$mapping_set_id);
    #            }
    #     }
  }
  elsif ($status == SAME_MAPPING){
    #seq_region_mapping has not change, nothing to do for the moment....

  }
  elsif ($status == NEW_MAPPING){
    $sth_mapping_set->execute($mapping_set_id,$schema_build);
    #there has been a seq_region change between releases, add a new mapping_set and the relation old_seq_region_id->new_seq_region_id
    my $previous_dbname = &get_previous_dbname($old_dbh,$db_name->[0],$release);
    $current_seq_region =  &read_seq_region($dbh,$db_name->[0]);
    $old_seq_region = &read_seq_region($old_dbh,$previous_dbname);
    #update the seq_region_mapping table with the old->new seq_region_id relation
    my $count = 0;
    foreach my $seq_region_name (keys %{$old_seq_region}){
      next if (!defined $current_seq_region->{$seq_region_name}); #the seq_region might have disappeared
      foreach my $coord_system_id (keys %{$old_seq_region->{$seq_region_name}}){
	next if (!defined $current_seq_region->{$seq_region_name}->{$coord_system_id}); #the coord_system might have been removed in current database
	next if ($old_seq_region->{$seq_region_name}->{$coord_system_id} == $current_seq_region->{$seq_region_name}->{$coord_system_id}); # if no change no need to write out
	$sth_seq_mapping->execute($old_seq_region->{$seq_region_name}->{$coord_system_id},$current_seq_region->{$seq_region_name}->{$coord_system_id},$mapping_set_id);
	$count++;
      }
    }
    print STDERR "Added $count seq_region_mapping entry\n\n";
  }
  else{
    throw("Mapping status not recognized by script: $status \n\n");
  }
}

#will for a given database, will return the seq_region_name->seq_region_id relation
sub read_seq_region{
    my $dbh = shift;
    my $dbname = shift;
    my %seq_region_hash;
    my $seq_region_id;
    my $seq_region_name;
    my $coord_system_id;
    my $sth = $dbh->prepare("SELECT seq_region_id, name, coord_system_id FROM $dbname.seq_region");
    $sth->execute();
    $sth->bind_col(1,\$seq_region_id);
    $sth->bind_col(2,\$seq_region_name);
    $sth->bind_col(3,\$coord_system_id);
    while ($sth->fetch){
	#there might be more than one assembly in the core database, thus we need the coord_system_id to remove ambiguity
	$seq_region_hash{$seq_region_name}{$coord_system_id} = $seq_region_id;
    }
    return \%seq_region_hash;
}

#method to check the status of the current core database: INITIAL_MAPPING, SAME_MAPPING and NEW_MAPPING are the possible states
sub mapping_status{
  my $dbh = shift;
  my $old_dbh = shift;
  my $dbname = shift;
  my $mapping_set_id_ref = shift;
  my $release = shift;
  
  #    my $sth_max_mapping = $dbh->prepare("select max(mapping_set_id) from $dbname.mapping_set");
  #    $sth_max_mapping->execute();
  #    ( $$mapping_set_id_ref ) = $sth_max_mapping->fetchrow_array();
  #    if (! $$mapping_set_id_ref){
  #	#the table is empty, first mapping
  #	$$mapping_set_id_ref = 1;
  #	return INITIAL_MAPPING;
  #    }
  #    else{
  #there is information, find out if it is the same mapping as previous release
  
  my $previous_dbname = &get_previous_dbname($old_dbh,$dbname,$release);
  if(!defined($previous_dbname)){
    print "No previous database present for $dbname so cannot do diff so will initialise with this as the first version of the database\n";
    $$mapping_set_id_ref = 1;
    return INITIAL_MAPPING;
  }
  my $cur_seq_region_size = &get_seq_region_size($dbh,$dbname);
  my $previous_seq_region_size = &get_seq_region_size($old_dbh,$previous_dbname);
  if ($cur_seq_region_size == $previous_seq_region_size){
    #if both tables have same size, SAME_MAPPING
    return SAME_MAPPING;
  }
  else{
    #if tables have different size, NEW_MAPPING
    $$mapping_set_id_ref++;
    return NEW_MAPPING;
  }	
  #}
}

#for a given database, returns the size of the seq_region_table
sub get_seq_region_size{
    my $dbh = shift;
    my $dbname = shift;
    my $sth_status = $dbh->prepare("show table status from $dbname like 'seq_region'") ;
    $sth_status->execute();
    my @table_status = $sth_status->fetchrow_array();
    return $table_status[6]; #return the size of the table
}

#will return the max mapping_set_id being used in the mapping_set table
sub get_max_mapping_set_id{
    my $dbh = shift;
    my $dbname = shift;

    my $sth_mapping = $dbh->prepare("select max(mapping_set_id) from mapping_set");
    $sth_mapping->execute();
    my ($max_mapping_set_id) = $sth_mapping->fetchrow_array();
    return $max_mapping_set_id;
}

#this method will return the name of the previous database to release for same species (assuming is present)
sub get_previous_dbname{
    my $dbh = shift;
    my $dbname = shift;
    my $release = shift;

    $dbname =~ /(^[a-z]+_[a-z]+_core_)/;

    my $previous_release_name = $1 . (--$release);
    my $previous_sth = $dbh->prepare("show databases like \'%$previous_release_name%\'");
    $previous_sth->execute();
    my ($previous_dbname) = $previous_sth->fetchrow_array();
    return $previous_dbname;
    
}

#for a standard ensembl database name, returns the release number and assembly
sub get_schema_and_build{
  my ($dbname) = @_;
  my @dbname = split/_/, $dbname;
  return join "_",$dbname[($#dbname -1)], $dbname[($#dbname)];
}

#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#The script populates a stable_id lookup database with all stable ids found in databases on a specified server for
#a specified db release.
#The stable ids are copied for objects listed in hash %group_objects

use strict;
use warnings;
use DBI qw( :sql_types );
use Getopt::Long;
use Data::Dumper;

my ( $rdbname, $version, $host, $port, $user, $pass );
my ( $lhost, $lport, $luser, $lpass, $test );

my $create       = 0;
my $create_index = 0;

my %group_objects = (
  core => {
    Exon        => 1,
    Gene        => 1,
    Transcript  => 1,
    Translation => 1,
    Operon      => 1,

# OperonTranscript => 1, # these are in transcript table anyway
    GeneArchive        => 1,
    TranscriptArchive  => 1,
    TranslationArchive => 1,
  },

# otherfeatures can be skipped for now as plant reuse the same ids in core and otherfeatures
  otherfeatures => {
    Gene        => 1,
    Transcript  => 1,
    Translation => 1,
  },

# not sure these are supposed to work - genetrees and families can't have species ids ..
# in fact in eg_23_76 there are no compara ids in the stable_id_lookup
  compara => {
    GeneTree => 1,
      Family   => 1,
  }
);

GetOptions(
  "lhost|lh=s"    => \$lhost,
  "lport=i"       => \$lport,
  "luser|lu=s"    => \$luser,
  "lpass|lp=s"    => \$lpass,
  "host|h=s"      => \$host,
  "port|p=i"      => \$port,
  "user|u=s"      => \$user,
  "pass=s"        => \$pass,
  "dbname|db=s"   => \$rdbname,
  "create!"       => \$create,
  "create_index!" => \$create_index,
  "version|v=s"   => \$version,
  "test!"         => \$test,
  "help", \&usage,

);

usage() if ( !defined $lhost || !defined $luser || !defined $lpass );
usage() if ( !defined $host  || !defined $port || !defined $user || !defined $version );

$rdbname ||= "ensembl_stable_ids_$version";
#$rdbname ||= "ensemblgenomes_stable_ids_$version";

my ( $readDB, $writeDB ) = init_db();

print "readDB credentials:\n";
print Dumper($readDB);

print "writeDB credentials:\n";
print Dumper($writeDB);

sub init_db {

# This is where the new stable_id database will be created
  my $writeDB_credentials = {
    user => defined($luser) ? $luser : $user,
    pass => defined($lpass) ? $lpass : $pass,
    host => defined($lhost) ? $lhost : $host,
    port => defined($lport) ? $lport : $port,
    dbname => $rdbname

  };

# This is where all the object dbs reside and queried against for stable ids
  my $readDB_credentials = {
    user   => $user,
    pass   => $pass,
    host   => $host,
    port   => $port,
    dbname => $rdbname
  };

  return ( $readDB_credentials, $writeDB_credentials );
}


create_db($writeDB) if $create;

my ( $dba_species, $lastSpeciesID ) = get_loaded_species($writeDB);
print("lastSpeciesID  + $lastSpeciesID\n");


my $new_species = {};

process_dbs($readDB);

create_index($writeDB) if $create_index;

sub create_index {
  my ($writeDB) = @_;

  my $dbname = $writeDB->{'dbname'};
  print "Creating index for $dbname\n";

  my $host = $writeDB->{'host'};
  my $port = $writeDB->{'port'};
  my $user = $writeDB->{'user'};
  my $pass = $writeDB->{'pass'};
  
  my $startAt = time;
  eval {
    my $cmd = "mysql -h $host";
    if ($port) {
      $cmd .= " -P $port";
    }
    $cmd .= " -u $user --password=$pass $dbname < ./sql/indices.sql";
    system($cmd) == 0 or die("error encountered when creating index for database $dbname\n");
  };

  if ($@) {
    die("An SQL error occured while creating database $dbname:\n$@");
  }
  my $took = time - $startAt;
  my $s    = $took % 60;
  my $m    = ( $took / 60 ) % 60;
  my $h    = $took / 3600;

  warn sprintf( "Indexed in %02d:%02d:%02d\n", $h, $m, $s );

}

sub create_db {
  my ($writeDB) = @_;

  my $dbname = $writeDB->{'dbname'};
  my $dbh = db_connect( 'test', $writeDB );
  print "Creating database $dbname\n";

  my $host = $writeDB->{'host'};
  my $port = $writeDB->{'port'};
  my $user = $writeDB->{'user'};
  my $pass = $writeDB->{'pass'};
  
  eval {

  $dbh->do("drop database if exists $dbname");
  $dbh->do("create database $dbname");

  my $cmd = "mysql -h $host";
  if ($port) {
    $cmd .= " -P $port";
  }
  $cmd .= " -u $user --password=$pass $dbname < ./sql/tables.sql";
  system($cmd) == 0 or die("error encountered when creating schema for database $dbname\n");

  $dbh->do("use $dbname");

  $dbh->do(
"INSERT INTO meta(species_id,meta_key,meta_value) VALUES (NULL,'schema_version','$version')"
  );

  };

  if ($@) {
    die("An SQL error occured while creating database $dbname:\n$@");
  }

  $dbh->disconnect();
}

sub process_dbs {
  my ($connectDB) = @_;

  my $host = $connectDB->{'host'};
  my $user = $connectDB->{'user'};
  my $pass = $connectDB->{'pass'};
  my $port = $connectDB->{'port'};

  my $counter_core  = 1;
  my $counter_other = 1;

  my $out;
  if ( defined $pass ) {
    $out = `mysql -h $host -P $port -u $user -p$pass -e 'show databases like "%$version%"'`;
  }
  else {
    $out = `mysql -h $host -P $port -u $user -e 'show databases like "%$version%"'`;
  }
  my @dbs = split /\n/, $out;

  my $startAt = time;

  foreach my $db (@dbs) {

    if ( $db =~ /([\w\_]+)_(core|otherfeatures)_([\d\_\w]+)/ )
    {
      my ( $species, $dbtype, $dbversion ) = ( $1, $2, $3 );
      print "SPECIES: $species\tDBTYPE $dbtype\tDBVERSION: $dbversion\n";
        if ( $dbversion =~ /^$version/ ) {
          if ( exists $dba_species->{$species} ) {
            if ( exists $dba_species->{$species}->{ID} ) {
              warn "* $species : LOADED\n";
              next;
            }
          }

        my $speciesOffset = $lastSpeciesID;
        if ( exists $new_species->{$species} ) {
          $speciesOffset = $new_species->{$species} - 1;    # this is the offset.
        }

        if ( $species =~ /_collection/ ) {
          add_collection_db( $db, $speciesOffset );
        }
        else {
          #print "\t\tGoing to add species $db   offset $speciesOffset\n";
          add_species_db( $db, $speciesOffset );
        }
        }
    }
    elsif ( $db =~ /([\w\_]+)_(compare)_([\d\_\w]+)/ ) {
      my ( $division, $dbtype, $dbversion ) = ( $1, $2, $3 );
      if ( $dbversion =~ /^$version/ ) {
        #add_compara_db($db);
      }
    }

  }

  my $took = time - $startAt;
  my $s    = $took % 60;
  my $m    = ( $took / 60 ) % 60;
  my $h    = $took / 3600;

  #warn sprintf( "Loaded in %02d:%02d:%02d\n", $h, $m, $s );
}

sub add_compara_db {
  my ($dbname) = @_;

# not sure what is supposed to happen to compara IDs - they do not have species_id
  warn "- Skipping $dbname\n";

}

sub add_species_db {
  my ( $dbname, $speciesOffset ) = @_;

  # 1 comes from species_id = 1 in meta table
  #print "From add_species_db\n";
  #warn "- Adding species $dbname (Species ID: ", 1 + $speciesOffset, ")\n";

  if ( $dbname =~ /([\w\_]+)_(core|otherfeatures)_([\d\_\w]+)/ ) {
    my ( $species, $dbtype, $dbversion ) = ( $1, $2, $3 );
    my $t1 = time;

    my $dba_read = db_connect( $dbname, $readDB );
    my $dba_write = db_connect( $writeDB->{"dbname"}, $writeDB );

    load_ids( $dba_read, $dba_write, $dbtype, $speciesOffset );

    if ( not exists $new_species->{$species} ) {
      load_species( $dba_read, $dba_write, $speciesOffset, $dbname );
      $new_species->{$species} = $lastSpeciesID;
    }

  $dba_read->disconnect();
  $dba_write->disconnect();

  #warn "+ Loaded in ", time - $t1, "s\n";
  }
}

sub add_collection_db {
  my ( $dbname, $speciesOffset ) = @_;

  warn "- Adding collection $dbname (from Species ID ", 1 + $speciesOffset,")\n";

  if ( $dbname =~ /([\w\_]+)_(core|otherfeatures)_([\d\_\w]+)/ ) {
    my ( $species, $dbtype, $dbversion ) = ( $1, $2, $3 );

    my $t1 = time;

    my $dba_read = db_connect( $dbname, $readDB );
    my $dba_write = db_connect( $writeDB->{"dbname"}, $writeDB );
    
    load_ids( $dba_read, $dba_write, $dbtype, $speciesOffset );

    if ( not exists $new_species->{$species} ) {
      load_species( $dba_read, $dba_write, $speciesOffset, $dbname );
      $new_species->{$species} = $lastSpeciesID;
    }

  $dba_read->disconnect();
  $dba_write->disconnect();

  #warn "+ Loaded in ", time - $t1, "s\n";
  }
}

sub load_species {
  my $dbh_read      = shift;
  my $dbh_write     = shift;
  my $speciesOffset = shift;
  my $dbname        = shift;

  my $sqlName = qq{SELECT species_id + $speciesOffset, meta_value FROM meta WHERE meta_key = "species.production_name"};
  my $sqlTaxon = qq{SELECT species_id + $speciesOffset, meta_value FROM meta WHERE meta_key = "species.taxonomy_id"};

  my $shash = {};

  my $sthN = $dbh_read->prepare($sqlName);
  $sthN->execute();
  while ( my ( $sid, $name ) = $sthN->fetchrow_array() ) {
    $shash->{$sid}->{Name} = $name;
      if ( $lastSpeciesID < $sid ) {
        $lastSpeciesID = $sid;
  }
  }
  $sthN->finish();

  my $sthT = $dbh_read->prepare($sqlTaxon);
  $sthT->execute();
  while ( my ( $sid, $taxid ) = $sthT->fetchrow_array() ) {
    $shash->{$sid}->{TaxID} = $taxid;
  }
  $sthT->finish();

  my @slist;
  my $insertSQL = qq{ INSERT INTO $writeDB->{'dbname'}.species (species_id, name, taxonomy_id) VALUES };
  my @tuples;

  foreach my $sid ( sort keys %{ $shash || {} } ) {
    push(@tuples, sprintf( q{(%s, %s, %s)}, $sid, $dbh_write->quote( $shash->{$sid}->{Name} ), $dbh_write->quote( $shash->{$sid}->{TaxID} ) )
  );
  }

  # Add the collection as well so if restart the script it does not load this collection again
  if ( $dbname =~ /([\w\_]+_collection)_(core|otherfeatures)_([\d\_\w]+)/ ) {
    my ( $species, $t, $v ) = ( $1, $2, $3 );
    push @tuples,
    sprintf( q{(%s, %s, 0)},
    $lastSpeciesID + 1,
    $dbh_write->quote($species) );
    $lastSpeciesID++;
  }

eval {
  $dbh_write->do( $insertSQL . join( ',', @tuples ) );
  if ($DBI::err) {
    warn $insertSQL . join( ',', @tuples ), "\n";
    die "ERROR: ", $DBI::errstr;
  }
  };
}

sub load_ids {

  my ( $dbh_read, $dbh_write, $dbtype, $speciesOffset ) = @_;
  #print "DbType : $dbtype     speciesOffset : $speciesOffset\n";

  my @stable_id_objects = keys %{ $group_objects{$dbtype} || {} };

  foreach my $object_name (@stable_id_objects) {
    my $object = lc($object_name);

    my $select_sql;

    if ( $object_name =~ /([A-Za-z]+)Archive/ ) {
      my $object = lc($1);

      #Note: Archive is not needed by ensemblgenomes as the stable_id_event is not populated
      #It is needed by ensembl. For the time being, we have to assume that the coord_system holds only one species with id of 1
      #In future ensembl also might need to support multiple databases
      my $species_id = 1;
      $select_sql =
      "SELECT DISTINCT old_stable_id, $species_id + $speciesOffset, '$dbtype', '$object' \
       FROM stable_id_event
       WHERE old_stable_id IS NOT NULL
       AND type = '$object'
       AND old_stable_id NOT IN (SELECT stable_id FROM $object)";
      $select_sql = $test ? $select_sql . " limit 10" : $select_sql;

      my $is_archive = 1;

      my $rows_inserted = build_insert_sql( $select_sql, $dbh_read, $dbh_write, $is_archive );
    }
    elsif ( $object_name =~ /Translation/ ) {
      my $sth = $dbh_read->prepare("SELECT COUNT(*) FROM $object");
      $sth->execute();
      my ($count) = $sth->fetchrow_array;
      if ($count) {
        if ( $dbtype eq 'core' ) {
          $select_sql =
          "SELECT DISTINCT o.stable_id, cs.species_id + $speciesOffset, '$dbtype', '$object_name' \
           FROM $object o \
           LEFT JOIN transcript t USING (transcript_id) \
           LEFT JOIN seq_region sr USING(seq_region_id) \
           LEFT JOIN coord_system cs USING(coord_system_id) \
           WHERE o.stable_id IS NOT NULL";

          $select_sql =  $test ? $select_sql . " limit 10" : $select_sql;

          my $rows_inserted = build_insert_sql( $select_sql, $dbh_read, $dbh_write );
        }
       elsif ( $dbtype eq 'otherfeatures' ) {

         $select_sql =
         "SELECT DISTINCT tl.stable_id, cs.species_id + $speciesOffset, '$dbtype', '$object_name' \
          FROM translation tl \
          LEFT JOIN transcript t USING (transcript_id) \
          LEFT JOIN analysis a USING (analysis_id) \
          LEFT JOIN seq_region sr USING(seq_region_id) \
          LEFT JOIN coord_system cs USING(coord_system_id) \
          WHERE logic_name like 'RefSeq_%' OR logic_name like 'CCDS_%'";

         $select_sql = $test ? $select_sql . " limit 10" : $select_sql;

         my $rows_inserted = build_insert_sql( $select_sql, $dbh_read, $dbh_write );
       }
      }
    }
    else {
      my $sth = $dbh_read->prepare("SELECT COUNT(*) FROM $object");
      $sth->execute();
      my ($count) = $sth->fetchrow_array;
      if ($count) {
        if ( $dbtype eq 'core' ) {
          $select_sql =
            "SELECT DISTINCT o.stable_id, cs.species_id + $speciesOffset, '$dbtype', '$object_name' \
             FROM $object o \
             LEFT JOIN seq_region sr USING(seq_region_id) \
             LEFT JOIN coord_system cs USING(coord_system_id) \
             WHERE o.stable_id is not NULL";

          $select_sql = $test ? $select_sql . " limit 10" : $select_sql;

          my $rows_inserted =
          build_insert_sql( $select_sql, $dbh_read, $dbh_write );
        }
        elsif ( $dbtype eq 'otherfeatures' ) {

          $select_sql =
            "SELECT DISTINCT o.stable_id, cs.species_id + $speciesOffset, '$dbtype', '$object_name' \
             FROM $object o \ 	
             LEFT JOIN analysis a USING (analysis_id) \
             LEFT JOIN seq_region sr USING(seq_region_id) \
             LEFT JOIN coord_system cs USING(coord_system_id) \
             WHERE logic_name like 'RefSeq_%' OR logic_name like 'CCDS_%'";

          $select_sql = $test ? $select_sql . " limit 10" : $select_sql;

          my $rows_inserted = build_insert_sql( $select_sql, $dbh_read, $dbh_write );
          }
        }
    }
  }

}

sub build_insert_sql {
  my ( $select_sql, $dbh_read, $dbh_write, $is_archive ) = @_;

  my $insert_sql = qq{INSERT INTO $writeDB->{'dbname'}.stable_id_lookup(stable_id, species_id, db_type, object_type) VALUES };

  if ($is_archive) {
    $insert_sql = qq{INSERT INTO $writeDB->{'dbname'}.archive_id_lookup(archive_id, species_id, db_type, object_type) VALUES };
  }

  my $import_sth = $dbh_read->prepare($select_sql);
  $import_sth->execute;
  my $rows = $import_sth->rows;

  my $start = 0;
  my $insert_values;
  my $max_rows_in_insert = 99999;    #Batch size of 100,000
  my $rows_in_insert     = 0;

  my @insert_container = ();
  my $while            = 0;
  while ( my $ref = $import_sth->fetchrow_arrayref ) {
    if ( $rows_in_insert < $max_rows_in_insert ) {
      $insert_values .= ',' if $start++;
      $insert_values .= '(' . ( join( ",", map { $dbh_read->quote($_) } @{$ref} ) ) . ')';
      $rows_in_insert++;
    }
    else {
      $insert_values .= ',' if $start++;
      $insert_values .= '(' . ( join( ",", map { $dbh_read->quote($_) } @{$ref} ) ) . ')';
      my $insert_sql_with_values = $insert_sql . $insert_values;
      push( @insert_container, $insert_sql_with_values );

      $insert_values  = undef;
      $rows_in_insert = 0;
      $start          = 0;
  }
}

  if ($insert_values) {
    my $insert_sql_with_values = $insert_sql . $insert_values;
    push( @insert_container, $insert_sql_with_values );
  }

  my $rows_inserted = 0;
  foreach my $insert_stmt (@insert_container) {
    my $affected = $dbh_write->do($insert_stmt);
    $rows_inserted += $affected;
  }

  #print "Number of rows fetched $rows \n Number of rows inserted into $writeDB->{'dbname'}.stable_id_lookup ", $rows_inserted, "\n";
  return $rows_inserted;

}

sub get_species {
  my $dbh    = shift;
  my $offset = shift;
  my $sql = qq{SELECT species_id, meta_value FROM meta WHERE meta_key = "species.production_name"};
  my $sth = $dbh->prepare($sql);
  $sth->execute();
  while ( my ( $sid, $name ) = $sth->fetchrow_array() ) {
    $dba_species->{$name}->{ID} = $sid + $offset;
  }
  $sth->finish();
}

sub db_connect {
  my ( $dbname, $connectDB ) = @_;

  $host = $connectDB->{'host'};
  $user = $connectDB->{'user'};
  $pass = $connectDB->{'pass'};
  $port = $connectDB->{'port'};

    my $dsn = "DBI:mysql:host=$host;";
    if ($port) {
      $dsn .= "port=$port;";
    }
    $dsn .= "database=$dbname";

    my $dbh = DBI->connect( $dsn, $user, $pass,
    { 'PrintError' => 1, 'RaiseError' => 1 } );

    if ( !$dbh ) {
      die "ERROR: $DBI::errstr";
    }

  return $dbh;
}

sub get_loaded_species {
  my $readDB = shift;

  my $shash = {};
  my $ssid  = 0;

  my $dbh = db_connect( $readDB->{'dbname'}, $readDB );

  my $sth = $dbh->prepare("SELECT species_id, name FROM species");
  $sth->execute();
  while ( my ( $sid, $name ) = $sth->fetchrow_array() ) {
    $shash->{$name}->{ID} = $sid;
    if ( $sid > $ssid ) {
      $ssid = $sid;
    }
  }
  $sth->finish();
  $dbh->disconnect();
  return ( $shash, $ssid );
}

sub usage {
  my $indent = ' ' x length($0);
  print <<EOF; exit(0);

The script populates a stable_id lookup database with all stable ids found in databases 
on a specified server for a specified db release.
Stable ids are copied for objects listed in hash %group_objects

The script is used for both ensembl and ensembl genomes

Options -lhost -lport -luser -lpass -version are mandatory and specify the credentials for the server on which a stable id lookup database exists or is to be created (if using option -create). 
If an argument for option -ldbname is not provided, the default name for the database wil be used: 'ensemblgenomes_stable_id_lookup_xx', where xx is the database release (option -version).

Options -host -user -port specify the credentials of the server(s) where stable ids are to be copied from.

To run the script cd into the directory where the script lives eg:
cd ensembl/misc-scripts/stable_id_lookup/


This command will create database ensembl_stable_ids_88 on server mysql-xxx-dev.ebi.ac.uk for release 88 databases found on mysql-ens-xxx.ebi.ac.uk:


eg: 
(with -test flag create a small subset database)
perl populate_stable_id_lookup.pl -lhost mysql-xxx-dev.ebi.ac.uk -luser xxxrw -lpass xxxx -lport 4484 -dbname ensembl_stable_ids_88 -create -host mysql-ens-xxx.ebi.ac.uk  -user ensro -port 4519 -version 88 -test

(without -test flag creates a full version)
perl populate_stable_id_lookup.pl -lhost mysql-xxx-dev.ebi.ac.uk -luser xxxrw -lpass xxxx -lport 4484 -dbname ensembl_stable_ids_88 -create -host mysql-ens-xxx.ebi.ac.uk  -user ensro -port 4519 -version 88


Usage:

  $0 -lhost host_name -lport port_number -luser user_name -lpass password -version db_version
  $indent [-create] [-dbname database_name] 
  $indent [-help]  
  
  -h|host              Database host where stable_ids are to be copied from (multiple hosts can be specified)

  -u|user              Database user where stable_ids are to be copied from (each host needs a user specified, 
                       if multiple -h|host options are given and fewer -u|user options are specified, 
                       the first user name will be used for the hosts where no user name was given)

  -port                Database port where stable_ids are to be copied from (if more than one host is specified 
                       multiple ports can be provided)

  -lh|lhost            Database host where stable_id lookup database exists or is to be created

  -lu|luser            Database user where stable_id lookup database exists or is to be created

  -lp|lpass            Database password where stable_id lookup database exists or is to be created

  -lport               Database port where stable_id lookup database exists or is to be created

  -v|version           EG version to match, e.g 88_77 OR Ensembl version to match 88

  -db|dbname           Database name for the stable id lookup database, e.g ensembl_stable_ids_88

  -create              Create the stable id lookup database using sql source ./sql/tables.sql

  -help                This message


EOF

}

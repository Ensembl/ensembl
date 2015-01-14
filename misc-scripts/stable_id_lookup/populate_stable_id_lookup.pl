#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;


my $lhost;
my $lport;
my $luser;
my $lpass;
my $ldbname;
my $create;
my $db_version;
my @host;
my @user;
my @port;

sub insert_stable_ids;
sub insert_species_id;

GetOptions( "lhost|lh=s" => \$lhost,
            "lport=i" => \$lport,
            "luser|lu=s" => \$luser,
            "lpass|lp=s" => \$lpass,
            "ldbname|ld=s" =>\$ldbname,
            "create!" => \$create,
            "db_version=i" => \$db_version,
            "host|h=s",\@host,
            "user|u=s",\@user,
            "port=i",\@port,
            "help" ,     \&usage,
        
);

usage() if (!defined $lhost || !defined $luser || !defined $lpass || !@host || !@user ); 


my $host_count = @host;
my $user_count = @user;
my $port_count = @port;

# if we have fewer user names specified than hosts copy user name from the first -u parameter
if ($user_count < $host_count) {
    for (my $i = $user_count; $i < $host_count; $i++) {
        push(@user,$user[0]);
    }
}

if ( (!@port) || ($port_count < $host_count) ) { 

    if (!defined $port[0]) {
        $port[0] = 3306;
    }

    for (my $i=1; $i<$host_count;$i++) {

        if (!defined $port[$i]) {
            push(@port,$port[0]);
        }
    } 
}

$db_version ||= software_version();
$ldbname ||= "ensembl_stable_ids_$db_version";
#$ldbname ||= "ensemblgenomes_stable_ids_$db_version";

my $registry = "Bio::EnsEMBL::Registry";


if ($host_count == 1) { 
    $registry->load_registry_from_db( -host => $host[0], -port => $port[0],-user => $user[0], -db_version => $db_version);
} else {

    my @server_array;

    for (my $i=0; $i < @host; $i++) {
        push @server_array, { -host => $host[$i], -user => $user[$i], -port => $port[$i]};
    }
        $registry->load_registry_from_multiple_dbs(@server_array); 

}
$registry->set_disconnect_when_inactive();

my @dbas = @{$registry->get_all_DBAdaptors()};

my $dbh; 
my $species_insert_sth;
my $species_sth;
my $stable_id_insert_sth;
my $archive_id_insert_sth;


if (@dbas) {

    #if any db adaptors exist (create and) connect to the stable id lookup database
    my $dsn = "DBI:mysql:host=$lhost;";
    if ($lport) {
        $dsn .= "port=$lport;";
    }
    if (!$create) {
        $dsn .= "database=$ldbname";
    }
    $dbh = DBI->connect( $dsn, $luser, $lpass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );
   
    if ($create) {
    print "Creating database $ldbname\n";

    eval {
        $dbh->do("drop database if exists $ldbname");
        $dbh->do("create database $ldbname");

        my $cmd = "mysql -h $lhost";
        if ($lport) {
            $cmd .= " -P $lport";
        }
        $cmd .= " -u $luser --password=$lpass $ldbname < ./sql/tables.sql";
        system($cmd) == 0 or die("error encountered when creating schema for database $ldbname\n");

        $dbh->do("use $ldbname");

        $dbh->do("INSERT INTO meta(species_id,meta_key,meta_value) VALUES (NULL,'schema_version',$db_version)");

    };

    if ($@) { 
        die("An SQL error occured while creating database $ldbname:\n$@");
    }


    }

    #statements used when populating the species table
    $species_sth = $dbh->prepare("SELECT species_id FROM species WHERE name = ?");
    $species_insert_sth = $dbh->prepare("INSERT INTO species(name,taxonomy_id) values (?,?)");
   
    #statements used when populating stable_id_lookup table
    $stable_id_insert_sth = $dbh->prepare("INSERT INTO stable_id_lookup VALUES(?,?,?,?)");
    $archive_id_insert_sth = $dbh->prepare("INSERT INTO archive_id_lookup VALUES(?,?,?,?)");
   
} else {
    die("No DBAdaptors found on ". join(',',@host) ." for db version $db_version\n");
}
 

my %group_objects = (
            core => {
                 Exon => 1,
                 Gene => 1,
                 Transcript => 1,
                 Translation => 1,
                 Operon => 1,
                 OperonTranscript => 1,
                 GeneArchive => 1,
                 TranscriptArchive => 1,
                 TranslationArchive => 1,
            },
            otherfeatures => {
                 Gene => 1,
                 Transcript => 1,
                 Translation => 1,
            },
                    );


#hash which stores species we have already processed
my %dba_species;


#populate stable_id_lookup table with stable ids 

while (my $dba = shift @dbas) {

    next if ( exists($dba_species{$dba->species()}{$dba->group()}) );

    my @stable_id_objects = keys %{$group_objects{$dba->group()}};

    my $species_id;

    if (@stable_id_objects) {

        my $species_name = $dba->species();
        $species_sth->bind_param( 1, $species_name, SQL_VARCHAR );
        $species_sth->execute();

        ($species_id) = $species_sth->fetchrow_array();

        if (!$species_id) {
            $species_id = insert_species_id($dba);
        }

        if ($species_id) {
            $dba_species{$dba->species()}{$dba->group()} = 1; 
        }     

    }

    foreach my $object_name (@stable_id_objects) {
    
        if ($object_name =~ /([A-Za-z]+)Archive/) {

            my $object = $1;
            my $lc_object = lc($object);

            my $dba_dbh = $dba->dbc->db_handle();
            my $archive_sql = qq(SELECT DISTINCT old_stable_id FROM stable_id_event
                                        WHERE old_stable_id IS NOT NULL
                                          AND type = '$lc_object'
                                          AND old_stable_id NOT IN (SELECT stable_id FROM $lc_object));
            my $archive_sth = $dba_dbh->prepare($archive_sql);
            $archive_sth->execute;
            my @archive_ids;

            while ((my $archive) = $archive_sth->fetchrow_array) {

               push (@archive_ids, $archive);

            }

            if (@archive_ids) {

                insert_ids(\@archive_ids, $dba, $object, $species_id, 'archive');

            }

        } elsif($dba->group eq 'otherfeatures') {

            my $import_sql;
            my $dba_dbh = $dba->dbc->db_handle();
            if ($object_name eq 'Translation') {

                $import_sql = qq(SELECT DISTINCT tl.stable_id FROM translation tl, transcript t, analysis a
                                             WHERE t.transcript_id = tl.transcript_id AND a.analysis_id = t.analysis_id
                                               AND (logic_name like 'RefSeq_%' OR logic_name like 'CCDS_%'));

            } else {

                my $object = lc($object_name);
                $import_sql = qq(SELECT DISTINCT stable_id FROM $object o, analysis a WHERE a.analysis_id = o.analysis_id AND (logic_name like 'RefSeq_%' OR logic_name like 'CCDS_%')); 

            }

            my $import_sth = $dba_dbh->prepare($import_sql);
            $import_sth->execute;
            my @import_ids;

            while ((my $import) = $import_sth->fetchrow_array) {

               push (@import_ids, $import);

            }

            if (@import_ids) {

                insert_ids(\@import_ids, $dba, $object_name, $species_id, 'stable');

            }

        } else {
            my $adaptor =  $dba->get_adaptor($object_name);
            my %stable_ids;
            my @archive_ids;

            if ($adaptor->can('list_stable_ids')) {
    
                %stable_ids = map { $_ => 1 } @{$adaptor->list_stable_ids()};
         
            } else {
    
                %stable_ids = map { ($_->stable_id() || '') => 1 } @{$adaptor->fetch_all()};
            }
    
    
            delete $stable_ids{''};
            my @stable_ids = keys %stable_ids;
     
            if (@stable_ids) {
    
                insert_ids(\@stable_ids, $dba, $object_name, $species_id, 'stable');
            }
        }
    }

}

$species_sth->finish() if ($species_sth);
$species_insert_sth->finish() if ($species_insert_sth);
$stable_id_insert_sth->finish() if ($stable_id_insert_sth);
$archive_id_insert_sth->finish() if ($archive_id_insert_sth);

$dbh->disconnect() if ($dbh);



sub insert_ids {

    my $ids = shift;
    my $dba = shift;
    my $object = shift;
    my $lookup_db_species_id = shift;
    my $id_type = shift;


    my @ids = @$ids;

    my @species_id;
    my @db_type;
    my @object_type;

    for (1..@ids) {
        push @species_id, $lookup_db_species_id;
        push @db_type, $dba->group();
        push @object_type, $object;
    }

    my $tuples;
    my @tuple_status;
    if ($id_type eq 'stable') {
      eval {
      $tuples = $stable_id_insert_sth->execute_array(
      { ArrayTupleStatus => \@tuple_status },
           \@ids,
           \@species_id,
           \@db_type,
           \@object_type,
       );
      };
    } elsif ($id_type eq 'archive') {
      eval {
      $tuples = $archive_id_insert_sth->execute_array(
      { ArrayTupleStatus => \@tuple_status },
           \@ids,
           \@species_id,
           \@db_type,
           \@object_type,
       );
      };
    }

    if ($tuples) {
        printf STDOUT "Successfully inserted %d %s ids for %s (%d), db type : %s, object type : %s\n", scalar @ids, $id_type, $dba->species(), $species_id[0], $dba->group(), $object;
    }
    else {
     for my $tuple (0..@ids-1) {
          my $status = $tuple_status[$tuple];
          $status = [0, "Skipped"] unless defined $status;
          next unless ref $status;
          printf STDERR "Failed to insert (%s, %s, %s, %s): %s\n",
          $ids[$tuple], $species_id[$tuple], $db_type[$tuple], $object_type[$tuple], $status->[1];
     }
    }

}


sub insert_species_id {

    my $dba = shift;
    my $species_name = $dba->species();

    #add species to the species table

    my $meta_container = $dba->get_adaptor('MetaContainer');
    my $taxonomy_id;
    if ($meta_container) {
       my $values = $meta_container->list_value_by_key('species.taxonomy_id'); 
       if ($values) {
           my @values = @$values;
           $taxonomy_id = $values[0];
       }
    }

    #add row to the species table
    $species_insert_sth->bind_param( 1, $species_name, SQL_VARCHAR );
    $species_insert_sth->bind_param( 2, $taxonomy_id, SQL_INTEGER );

    $species_insert_sth->execute();
    my $species_id = $dbh->last_insert_id( undef, undef, 'species', 'species_id' ); 

    if (!$species_id) {
        die("Failed to insert row for species $species_name\n");
    }
    return $species_id;

}


sub usage {
  my $indent = ' ' x length($0);
  print <<EOF; exit(0);

The script populates a stable_id lookup database with all stable ids found in databases 
on a specified server (or servers) for a specified db release.
Stable ids are copied for objects listed in hash %group_objects

Options -lhost -luser -lpass are mandatory and specify the credentials for the server on which a stable id lookup database exists or is to be created (if using option -create). If an argument for option -ldbname is not provided, the default name for the database wil be used: 'ensembl_stable_id_lookup_xx', where xx is the database release (option -db_version).

Options -host -user -port specify the credentials of the server(s) where stable ids are to be copied from.

To run the script cd into the directory where the script lives eg:
cd ensembl/misc-scripts/stable_id_lookup/


This command will create database ensembl_stable_id_lookup_67 on server ens-staging1 and will copy stable ids from databases for release 67 found on ens-staging1 and ens-staging2:

populate_stable_id_lookup.pl -lhost ens-staging1 -luser ensadmin -lpass xxxx -create -db_version 67 -host ens-staging1 -host ens-staging2 -user ensro


Usage:

  $0 -lhost host_name -luser user_name -lpass password 
  $indent [-ldbname database_name] [-lport port_number] 
  $indent -host host_name [-host host_name2] -user user_name [-user user_name2] 
  $indent [-port port_number [-port port_number2]]
  $indent [-create] [-db_version]
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

  -ld|ldbname          Database name for the stable id lookup database   

  -create              Create the stable id lookup database using sql source ./sql/tables.sql

  -db_version          If not specified, software_version() returned by the ApiVersion module will be used

  -help                This message


EOF

}

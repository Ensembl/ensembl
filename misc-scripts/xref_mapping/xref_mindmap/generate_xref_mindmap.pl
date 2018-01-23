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


#The script generates a FreePlane mind map file of xrefs for a given species.
#To view the generated file, download free mind mapping software Freeplane from http://freeplane.sourceforge.net

use strict;
use warnings;
use DBI qw( :sql_types );
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
use HTML::Entities;


use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin";
}

my $xhost;
my $xport;
my $xuser;
my $xpass;
my $db_version;
my $host;
my $user;
my $port;
my $species;
my $dontdrop;



GetOptions( "xhost|xh=s" => \$xhost,
            "xport=i" => \$xport,
            "xuser|xu=s" => \$xuser,
	    "xpass|xp=s" => \$xpass,
	    "host|h=s",\$host,
	    "user|u=s",\$user,
	    "port=i",\$port,
	    "species|s=s",\$species,
	    "db_version=i",\$db_version,
	    "dontdrop|d", \$dontdrop,
	    "help" ,     \&usage,
	    
);

usage() if (!defined $xhost || !defined $xuser || !defined $xpass || !defined $host || !defined $user || !defined $species); 



$port ||= 3306;
$xport ||= 3306;
$db_version ||= software_version();

my $registry = "Bio::EnsEMBL::Registry";


$registry->load_registry_from_db( -host => $host, -port => $port,-user => $user, -db_version => $db_version, -species => $species);


my $dba = $registry->get_DBAdaptor($species, "core");

if (!$dba) {
    die("$species database version $db_version not found on host $host port $port");
}

my $dbh; 

#create and connect to the xref_mindmap database
my $dsn = "DBI:mysql:host=$xhost;";
if ($xport) {
   $dsn .= "port=$xport;";
}  
my $xdbname = 'xref_mindmap_'.$species.'_'.$db_version;


$dbh = DBI->connect( $dsn, $xuser, $xpass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );


my @objects_with_xrefs = ('Gene', 'Transcript', 'Translation');
    
my %object_display_names = ('Gene' => 'Gene', 'Transcript' => 'Transcript', 'Translation' => 'Protein');


print STDOUT "Creating database $xdbname\n";

eval {
     $dbh->do("drop database if exists $xdbname");
     $dbh->do("create database $xdbname");

     my $cmd = "mysql -h $xhost";
     if ($xport) {
		$cmd .= " -P $xport";
     }
     $cmd .= " -u $xuser --password=$xpass $xdbname <  $SERVERROOT/sql/tables.sql";
     system($cmd) == 0 or die("error encountered when creating schema for database $xdbname\n");

     $dbh->do("use $xdbname");

};

if ($@) { 
     die("An SQL error occured while creating database $xdbname:\n$@");
}



#copy the external_db information

my $external_db_info_sth = $dba->dbc()->prepare("select external_db_id, db_name, db_display_name from external_db");
my $insert_external_db_info_sth = $dbh->prepare("insert into external_db_type (external_db_id, db_name, db_display_name) values (?,?,?)");


$external_db_info_sth->execute();

while ( my ($external_db_id, $db_name, $db_display_name) = $external_db_info_sth->fetchrow_array() ) {

    $insert_external_db_info_sth->execute($external_db_id, $db_name, $db_display_name);

}
$external_db_info_sth->finish();
$insert_external_db_info_sth->finish();

#populate db types in external_db_type

my $cmd = "mysql -h $xhost";
if ($xport) {
    $cmd .= " -P $xport";
}
$cmd .= " -u $xuser --password=$xpass $xdbname <  $SERVERROOT/sql/update.sql";
system($cmd) == 0 or die("error encountered when updating database $xdbname\n");

 
#get distinct db_names for an object type
my $external_db_sth = $dba->dbc()->prepare("select distinct external_db_id from object_xref join xref using(xref_id) join external_db using(external_db_id) where ensembl_object_type = ? and info_type != 'UNMAPPED' order by db_name");

#get xref types for an object type and external db
my $xref_types_sth = $dba->dbc()->prepare("select distinct info_type from object_xref join xref using(xref_id) where ensembl_object_type = ? and external_db_id = ? and info_type != 'UNMAPPED'");

#get master xref db_names
my $master_db_sth =  $dba->dbc()->prepare("select distinct e.external_db_id from object_xref ox join xref x1 on (ox.xref_id = x1.xref_id) join dependent_xref on (x1.xref_id = dependent_xref_id) join xref x2 on (master_xref_id = x2.xref_id) join external_db e on (x2.external_db_id = e.external_db_id) where ensembl_object_type = ? and  x1.external_db_id = ? and x1.info_type != 'UNMAPPED'");

#get species from which xrefs were projected
my $projected_from_species_sth = $dba->dbc()->prepare("select distinct substr(info_text,1,locate(' ', info_text, locate(' ', info_text)+1 )) from object_xref join xref using(xref_id) where info_text like 'from%' and ensembl_object_type = ? and external_db_id = ? and info_type = 'PROJECTION'");


#xref_mindmap db statements

my $insert_object_xref_linkage = $dbh->prepare("insert into object_xref_linkage values(?,?,?,?,?)");

my $insert_object_external_db_node = $dbh->prepare("insert into object_external_db_node values (?,?,?)");

my $link_type_sth = $dbh->prepare("select link_type_id from link_type where link_type = ?");


foreach my $object (@objects_with_xrefs) {

    $external_db_sth->execute($object);
    my @db_ids;
    while (my ($db_id) = $external_db_sth->fetchrow_array() ) {
	push @db_ids, $db_id;
    }

    foreach my $db_id (@db_ids) {
	$xref_types_sth->execute($object, $db_id);

	#store the xref db node
	$insert_object_external_db_node->execute($object,$db_id, $object.$db_id);

	while (my ($link_type) = $xref_types_sth->fetchrow_array() ) {

	    #get link_type_id
	    $link_type_sth->execute($link_type);
	    my ($link_type_id) = $link_type_sth->fetchrow_array();

	    if ($link_type_id) {

		if ($link_type eq 'DEPENDENT') {

		   $master_db_sth->execute($object,$db_id);
		   while (my ($master_db_id) = $master_db_sth->fetchrow_array() ) {

		       $insert_object_xref_linkage->execute($db_id,$object,$link_type_id, $master_db_id, undef);
		   }

		} elsif ($link_type eq 'PROJECTION') {
		
		    $projected_from_species_sth->execute($object, $db_id);
		    while (my ($linked_node_text) = $projected_from_species_sth->fetchrow_array()) {
			$insert_object_xref_linkage->execute($db_id,$object,$link_type_id, undef, $linked_node_text);
		    }

		} else {
		    
		    $insert_object_xref_linkage->execute($db_id,$object,$link_type_id, undef, undef);
		
		}
	    
	    }
	  
	}
	
    }

}


#populate 'generated_from' link types in object_xref_linkage - this is for xrefs to do with gene and transcript naming

my $naming_dbs_sth = $dbh->prepare("select external_db_id, db_name from external_db_type where db_name like ?");
my $db_id_sth = $dbh->prepare("select external_db_id from external_db_type where db_name = ?");

$link_type_sth->execute('GENERATED_FROM');
my ($link_type_id) = $link_type_sth->fetchrow_array();  

foreach my $object (@objects_with_xrefs) {

    $naming_dbs_sth->execute('%'.$object.'\_name');
    while (my ($db_id, $db_name) = $naming_dbs_sth->fetchrow_array() ) {

	if ($db_name =~ m/(.*)\_$object\_name/i) {

	    my $linked_db_name = $1;
	    if ($linked_db_name) {
	
	       #find the linked db id
	       $db_id_sth->execute($linked_db_name);
	       my ($linked_db_id) = $db_id_sth->fetchrow_array();

	       if ($linked_db_id) {
		    $insert_object_xref_linkage->execute($db_id,$object,$link_type_id, $linked_db_id, undef);
	       }
		
	    }
	}
    }
    
    foreach my $name ( ('Vega','Ensembl')) {
	$naming_dbs_sth->execute('Clone\_based\_'.$name.'\_'.$object);
	my ($db_id, $db_name) = $naming_dbs_sth->fetchrow_array();	
	if ($db_id) {
	    $insert_object_xref_linkage->execute($db_id,$object,$link_type_id, undef, "$name clone name");
	}
    }

}
$naming_dbs_sth->finish();


#check if there are any rows in the interpro table

my $interpro_sth = $dba->dbc()->prepare("select count(1) from interpro");

$interpro_sth->execute();
my ($interpro_count) = $interpro_sth->fetchrow_array();

if ($interpro_count > 0) {

    my $interpro = 'Interpro';

    $db_id_sth->execute($interpro);
    my ($interpro_db_id) = $db_id_sth->fetchrow_array();

    $link_type_sth->execute('PROTEIN_FEATURES');
    my ($link_type_id) = $link_type_sth->fetchrow_array();   

     my $object = 'Translation';
    $insert_object_external_db_node->execute($object,$interpro_db_id, $object.$interpro_db_id);

    $insert_object_xref_linkage->execute($interpro_db_id,$object,$link_type_id, undef, undef);
}

$interpro_sth->finish();
$db_id_sth->finish();

$external_db_sth->finish();
$xref_types_sth->finish();
$master_db_sth->finish();
$projected_from_species_sth->finish();

$insert_object_xref_linkage->finish();
$insert_object_external_db_node->finish();
$link_type_sth->finish();



#create the .mm file based on info from xref_mindmap db

my $file_name = $species .'_xrefs_'. $db_version .'.mm';

open( FH, ">$file_name" ) or die("Can't open $file_name\n");

print STDOUT "Writing to file $file_name\n";


my $header = <<HEADER;
<map version='0.9.0'>
<!--To view this file, download free mind mapping software Freeplane from http://freeplane.sourceforge.net -->
<node TEXT='External References' ID='ID_000000001' COLOR='#18898b' STYLE='fork'>
<font NAME='Liberation Sans' SIZE='12' BOLD='true'/>
<hook NAME='MapStyle' max_node_width='600'/>
<edge STYLE='bezier' COLOR='#808080' WIDTH='thin'/>
HEADER

print FH $header;


my $category_and_db_name_nodes_sth = $dbh->prepare("select db_type, external_db_id, db_display_name, mindmap_tag_id from object_external_db_node join external_db_type using(external_db_id) join db_type using(db_type_id) where ensembl_object_type = ? order by db_type, db_display_name;"); 

my $db_linkage_types_sth = $dbh->prepare("select link_description, linked_external_db_id, linked_node_text from object_xref_linkage join link_type using(link_type_id) where ensembl_object_type = ? and external_db_id = ? and link_description is not null order by link_description, linked_external_db_id, linked_node_text");

#this is used to link to dbs from the closest level (e.g. if Orphanet xrefs for a Gene are dependent on HGNC xrefs
#we will try to link to the Gene HGNC node then transcript etc..
 
my $object_distance_sth = $dbh->prepare("select to_object from object_distance where from_object = ? order by distance");

my $node_id_sth = $dbh->prepare("select mindmap_tag_id from object_external_db_node where ensembl_object_type = ? and external_db_id = ?");


my $node_id_count = 1;

foreach my $object (@objects_with_xrefs) {
    $node_id_count++;
    my $node_id = 'ID_'.$node_id_count;
    #write out the object node
    my $object_node = "<node TEXT='".$object_display_names{$object}."' POSITION='right' ID='".$node_id."' COLOR='#cc3300' STYLE='fork'>
<font NAME='Liberation Sans' SIZE='12' BOLD='true'/>
<edge STYLE='bezier' COLOR='#808080' WIDTH='thin'/>
";
    print FH $object_node;
  
    $category_and_db_name_nodes_sth->execute($object);

    my $last_type = 'none';
    while (my ($db_type, $db_id, $db_display_name, $tag_id) = $category_and_db_name_nodes_sth->fetchrow_array() ) {

	if ($db_type ne $last_type) {

	    if ($last_type ne 'none') {
		print FH "</node>\n";  #db type node end tag
	    }
	    $node_id_count++;
	    $node_id = 'ID_'.$node_id_count;
	    #write out new db type node
	    my $db_type_node = "<node TEXT='".encode_entities($db_type)."' ID='".$node_id."' COLOR='#669900'>
<font NAME='Liberation Sans' SIZE='12' BOLD='true'/>
";
	    print FH $db_type_node;
	}

	#write out the db name node
	$node_id_sth->execute($object,$db_id);
	($node_id) = $node_id_sth->fetchrow_array(); 
	my $db_name_node = "<node TEXT='".encode_entities($db_display_name)."' ID='".$node_id."'>
";
	print FH $db_name_node;

	$db_linkage_types_sth->execute($object,$db_id);

	my $last_link_desc = 'none';
	while (my ($link_description, $linked_db_id, $linked_node_text) 
	       = $db_linkage_types_sth->fetchrow_array() ){

	    if ($link_description ne $last_link_desc) {
		if ($last_link_desc ne 'none') {
		    print FH "</node>\n";  #link type node end tag
		}
		#write out the link node
		$node_id_count++;
		$node_id = 'ID_'.$node_id_count;
		my $link_node = "<node TEXT='".encode_entities($link_description)."' ID='".$node_id."'>
";
		print FH $link_node;

	    } 

	    if ($linked_db_id) {

		#find the node id to link to 
		my @closest_level;
		$object_distance_sth->execute($object);
	       
		while (my ($object_type) = $object_distance_sth->fetchrow_array() ) {
		    push @closest_level, $object_type;
		} 

		my $linked_node_id;
		foreach my $object_type (@closest_level) {
		    
		    $node_id_sth->execute($object_type,$linked_db_id);
		    ($linked_node_id) = $node_id_sth->fetchrow_array();
		    last if ($linked_node_id);
		}

		if ($linked_node_id) {
		    #write out arrow link to the node
		    my $arrow_link_node = "<arrowlink DESTINATION='".$linked_node_id."' STARTARROW='NONE' ENDARROW='DEFAULT'/>
";
		    print FH $arrow_link_node;
		}

	    } elsif ($linked_node_text) {

		#write out the node
		$node_id_count++;
		$node_id = 'ID_'.$node_id_count;
		my $node = "<node TEXT='".encode_entities($linked_node_text)."' ID='".$node_id."'/>
";
		print FH $node;		

	    }

	    $last_link_desc = $link_description;


	}
	if ($last_link_desc ne 'none') {
	    print FH "</node>\n";  #link type node end tag
	}

	
	print FH "</node>\n";  #db name node end tag

	$last_type = $db_type;
    }
    if ($last_type ne 'none') {
	print FH "</node>\n";  #db type node end tag
    }


    print FH "</node>\n" #object node end tag;

}

$category_and_db_name_nodes_sth->finish();
$db_linkage_types_sth->finish();
$object_distance_sth->finish();
$node_id_sth->finish();


#external references end tag and map end tag
my $footer = <<FOOTER;
</node>
</map>
FOOTER

print FH $footer;

close(FH);

if (!$dontdrop) {

    print STDOUT "Dropping database $xdbname\n";
    eval {
	$dbh->do("drop database if exists $xdbname");
    };
    if ($@) { 
	die("An SQL error occured while dropping database $xdbname:\n$@");
    }
}

$dbh->disconnect();




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

  $0 -xhost host_name -xuser user_name -xpass password 
  $indent [-xport port_number] 
  $indent -host host_name -user user_name
  $indent [-port port_number] [-db_version]
  $indent -species species_name
  $indent [-dontdrop] 
  $indent [-help]  

  

  -h|host              Database host where the species db lives

  -u|user              Database user where the species db lives

  -port                Database port the species db lives

  -xh|xhost            Database host where the xref mindmap db is to be created

  -xu|xuser            Database user where the xref mindmap db is to be created

  -xp|xpass            Database password where the xref mindmap db is to be created

  -xport               Database port where the xref mindmap db is to be created

  -s|species           Species name to generate xref mindmap for

  -d|dontdrop          Don\'t drop the xref_mindmap db at the end of the script

  -db_version          If not specified, software_version() returned by the ApiVersion module will be used

  -help                This message


EOF

}

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

# this script will repopulate the meta coord table
# it makes an attempt to pick up all the right tables
# (all those which have seq_region_id, seq_region_start and seq_region_end)

# normally the API will populate the table ...


use strict;
use warnings;
use DBI;

use Getopt::Long;

my ($host, $user, $pass, $port, $dbname); #ensembl core db
GetOptions('host=s'   => \$host,
	   'user=s'   => \$user,
	   'pass=s'   => \$pass,
	   'port=i'   => \$port,
	   'dbname=s' => \$dbname,
	  );


my $dsn = "DBI:mysql:host=$host;dbname=$dbname";
if( $port ) {
  $dsn .= ";port=$port";
}

my $db = DBI->connect( $dsn, $user, $pass );

$db->do( "delete from meta_coord" );

my $res = $db->selectall_arrayref( "show tables" );

my @tables = map { $_->[0] } @$res;

my %need_cols = ( "seq_region_id" => 1,
		  "seq_region_start" => 1,
		  "seq_region_end" => 1 );

for my $tablename ( @tables ) {
  # skip empty tables
  $res = $db->selectall_arrayref( "select count(*) from $tablename" );
  next if( $res->[0]->[0] == 0 );

  $res = $db->selectall_arrayref( "desc $tablename" );
  my @columns = map { $_->[0] } @$res;
  if( 3 == scalar( grep { exists $need_cols{$_}} @columns )) {
    meta_coord_query( $db, $tablename );
  }
}




sub meta_coord_query {
  my $db = shift;
  my $tablename = shift;

  $db->do( qq{
	      INSERT INTO meta_coord
	      SELECT '$tablename', sr.coord_system_id,
	      MAX( f.seq_region_end -  f.seq_region_start + 1)
	      FROM   $tablename f, seq_region sr
	      WHERE  sr.seq_region_id = f.seq_region_id
	      GROUP BY sr.coord_system_id
	      } );
}



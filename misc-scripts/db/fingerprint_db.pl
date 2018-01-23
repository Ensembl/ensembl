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

use strict;
use warnings;


use DBI;
use Getopt::Long;

my ($host, $port, $dbname, $user, $pass, @tables, $limit);

GetOptions('user=s'       => \$user,
	   'pass=s'       => \$pass,
	   'host=s'       => \$host,
	   'port=i'       => \$port,
	   'dbname=s'     => \$dbname,
	   'tables=s'     => \@tables,
	   'limit=s'      => \$limit,
	   'help'         => sub { usage(); exit(0); });

@tables = split(/,/,join(',',@tables));

if (!$user || !$host || !$dbname) {

  usage();
  exit(1);

}

my $dbi = DBI->connect( "DBI:mysql:host=$host:port=$port;database=$dbname", $user, $pass,
			{'RaiseError' => 1}) || die "Can't connect to database\n";

# use all if not specified
if (!@tables) {
  my $sth = $dbi->prepare("SHOW TABLES;");
  $sth->execute();
  foreach my $row (@{$sth->fetchall_arrayref()}) {
    push @tables, $row->[0];
  }
}

print "Table\tDescription checksum\tData checksum\n";

foreach my $table (@tables) {

  print $table . "\t";

  my $passopt = $pass ? "-p$pass" : "";

  my $key_sth = $dbi->prepare("SHOW INDEX FROM $table");
  $key_sth->execute();
  my ($primary_key, $first_col);
  while (my @cols = $key_sth->fetchrow()) {
    $first_col = $cols[4];
    if ($cols[2] eq "PRIMARY") {
      $primary_key = $cols[4];
      last;
    }
  }

  my $key = $primary_key;;
  if (!$key) {
    # print "Can't get primary key for table $table, using $first_col as key\n";
    $key = $first_col;
  }

  # description - note use of tr -d to strip newline
  system ("mysql -u $user -h $host $passopt -P $port -e \'DESC $table\' $dbname | md5sum | tr -d '\n'");

  print "\t";

  # data
  my $limitclause = $limit ? " LIMIT $limit" : "";
  my $cmd = "\'SELECT * FROM $table ORDER BY $key $limitclause;desc $table;\'";
  system ("mysql -u $user -h $host $passopt -P $port -e $cmd $dbname | md5sum");

}



sub usage {

  print << "EOF";

  fingerprint_db.pl -user {user} -pass {password} -host {host} -port {port} -dbname {database} -tables {tables} -limit {rows}

  Generates an MD5 checksum for the description and data of all or some tables of a database.

  Argument to tables may be comma-separated, or multiple -tables arguments may be used.
  If no -tables argument is given all tables are analysed.

  Use -limit to limit the number of rows the comparison is based on; less accurate but faster for large tables.

  A single fingerprint for the whole database can be computed by piping the output of this script into md5sum.

EOF

}

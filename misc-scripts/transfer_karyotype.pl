#!/usr/bin/env perl
#
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

# Transfer karyotype data from old -> new database, scaling appropriately.

use strict;
use warnings;


use DBI;
use Getopt::Long;

my ( $oldhost, $olduser, $oldpass, $oldport, $olddbname, $newhost, $newuser, $newpass, $newport, $newdbname  );


GetOptions( "oldhost=s",   \$oldhost,
	    "olduser=s",   \$olduser,
	    "oldpass=s",   \$oldpass,
	    "oldport=i",   \$oldport,
	    "olddbname=s", \$olddbname,
	    "newhost=s",   \$newhost,
	    "newuser=s",   \$newuser,
	    "newpass=s",   \$newpass,
	    "newport=i",   \$newport,
	    "newdbname=s", \$newdbname);

usage() if (!$oldhost);

my $old_db = DBI->connect("DBI:mysql:host=$oldhost;dbname=$olddbname;port=$oldport", $olduser, $oldpass);

my $new_db = DBI->connect("DBI:mysql:host=$newhost;dbname=$newdbname;port=$newport", $newuser, $newpass);

# check for rows in new table
my $chk_sth = $new_db->prepare("SELECT COUNT(*) FROM karyotype");
$chk_sth->execute();
my $count = ($chk_sth->fetchrow_array())[0];
if ($count > 0) {
  print STDERR "Karyotype table in $newdbname should be empty but has $count rows - delete and re-run\n";
  exit(1);
}

my $old_sth = $old_db->prepare("SELECT sr.name, cs.name, sr.length, k.seq_region_start, k.seq_region_end, k.band, k.stain FROM seq_region sr, coord_system cs, karyotype k WHERE sr.coord_system_id=cs.coord_system_id AND sr.seq_region_id=k.seq_region_id");
$old_sth->execute();

my ($old_sr_name, $old_cs_name, $old_sr_length, $old_k_start, $old_k_end, $band, $stain);
$old_sth->bind_columns(\$old_sr_name, \$old_cs_name, \$old_sr_length, \$old_k_start, \$old_k_end, \$band, \$stain);

my $new_sth = $new_db->prepare('SELECT sr.seq_region_id, sr.length FROM seq_region sr, coord_system cs WHERE sr.name=? and cs.name=? AND sr.coord_system_id=cs.coord_system_id AND cs.attrib like "%default_version%" ');
my ($new_sr_id, $new_sr_length);

my $insert_sth = $new_db->prepare("INSERT INTO karyotype (seq_region_id, seq_region_start, seq_region_end, band, stain) VALUES(?,?,?,?,?)");

my $count;

while ($old_sth->fetch()) {

  # get matching seq region from new database & calculate scaling factor
  $new_sth->execute($old_sr_name, $old_cs_name);
  $new_sth->bind_columns(\$new_sr_id, \$new_sr_length);
  my $scale_factor;
  if ($new_sth->fetch()) {

    $scale_factor = ($new_sr_length/$old_sr_length);
    my $new_k_start = int($old_k_start * $scale_factor) || 1;
    my $new_k_end = int($old_k_end * $scale_factor);
    if ($old_k_end == $old_sr_length) {
      $new_k_end = $new_sr_length;
    }

    # Add new entry to new karyotype table
    $insert_sth->execute($new_sr_id, $new_k_start, $new_k_end, $band, $stain) || die "Error inserting into new karyotype table";
    $count++;

  } else {
    warn("Can't get new seq_region ID corresponding to $old_cs_name:$old_sr_name\n");
  }

}

print "Inserted $count rows into $newdbname.karyotype\n";


sub usage {

  print<<EOF;
Transfer karyotype data from old to new database, scaling appropriately.

perl transfer_karyotype.pl -oldhost ... -olduser ... -oldpass ... -oldport ... -olddbname ... -newhost ... -newuser ... -newpass ... -newport ... -newdbname

EOF

  exit(1);

}

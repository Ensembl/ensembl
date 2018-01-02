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

# Compare xrefs between 2 databases by source

use strict;
use warnings;


use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Getopt::Long;

my ($old_host, $old_user, $old_pass, $old_port, $old_dbname, $new_host, $new_user, $new_pass, $new_port, $new_dbname, $source);

GetOptions( "old_host=s",   \$old_host,
	    "old_user=s",   \$old_user,
	    "old_pass=s",   \$old_pass,
	    "old_port=i",   \$old_port,
	    "old_dbname=s", \$old_dbname,
	    "new_host=s",   \$new_host,
	    "new_user=s",   \$new_user,
	    "new_pass=s",   \$new_pass,
	    "new_port=i",   \$new_port,
	    "new_dbname=s", \$new_dbname,
	    "source=s",     \$source,
	    "help",         \&usage
	  );

my $old_dbi = DBI->connect("dbi:mysql:host=$old_host;port=$old_port;database=$old_dbname",
			   $old_user, $old_pass,
			   {'RaiseError' => 1}) || die "Can't connect to database $old_dbname";

my $new_dbi = DBI->connect("dbi:mysql:host=$new_host;port=$new_port;database=$new_dbname",
			   $new_user, $new_pass,
			   {'RaiseError' => 1}) || die "Can't connect to database $new_dbname";


# if one source specified, use that, otherwise loop over them all

if ($source) {

  compare($source);

} else {

  my $old_extdb_sth = $old_dbi->prepare("SELECT db_name FROM external_db WHERE db_name NOT LIKE 'AFFY%' AND status LIKE 'KNOWN%'");
  my $ext_db;
  $old_extdb_sth->execute();
  $old_extdb_sth->bind_columns(\$ext_db);
  while ($old_extdb_sth->fetch()) {

    compare($ext_db);

  }

}



sub compare {

  my $source = shift;

  # Read & cache all old stable_id-display_xref mappings

  #print "Caching old gene stable ID - $source display_xref mappings from $old_dbname\n";

  my $sql = "SELECT g.stable_id, x.dbprimary_acc FROM xref x, gene g, external_db e WHERE g.display_xref_id=x.xref_id AND e.external_db_id=x.external_db_id AND e.db_name='" . $source . "'";

  #print "\n\n$sql\n\n";
  my $old_sth = $old_dbi->prepare($sql);

  my ($stable_id, $accession);
  my %old_stable_id_2_accession;
  $old_sth->execute();
  $old_sth->bind_columns(\$stable_id, \$accession);
  while ($old_sth->fetch()) {

    $old_stable_id_2_accession{$stable_id} = $accession;

  }

  my ($match, $mismatch, $added) = (0,0,0);

  # Compare these with new ones
  my $new_sth = $new_dbi->prepare($sql);

  $new_sth->execute();
  $new_sth->bind_columns(\$stable_id, \$accession);
  while ($new_sth->fetch()) {

    if (exists($old_stable_id_2_accession{$stable_id})) {
      if ($old_stable_id_2_accession{$stable_id} eq $accession) {
	$match++;
      } else {
	$mismatch++;
      }
    } else {
      $added++;
    }
  }


  printf "%30s\tmatch: %d\tmismatch: %d\tadded: %d\n", $source, $match, $mismatch, $added  if ($match || $mismatch || $added);

}



sub usage {

  print <<EOF;
Usage: perl compare_xrefs_by_source.pl

Compare assigned display xrefs between two databases, optionally by source.

Options:

       -old_host -old_user -old_pass -old_port -old_dbname  : old database connection details
       -new_host -new_user -new_pass -new_port -new_dbname  : new database connection details

       -source : source to compare (e.g. RefSeq_dna). If no source is specified, all sources
                 from old database are checked.

Note only non-Affy sources of status "KNOWN" or "KNOWNXREF" are compared, and if no differences
are detected, nothing is reported.

EOF

  exit(0);

}

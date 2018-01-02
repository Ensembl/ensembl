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

use Bio::EnsEMBL::DBSQL::DBConnection;
use Getopt::Long;
use Bio::EnsEMBL::Utils::CliHelper;

my $help = 0;
my ( $host, $port, $user, $pass, $dbpattern );

$port = '3306';

my @table_names = qw(
  assembly_exception
  density_feature
  ditag_feature
  dna_align_feature
  exon
  gene
  intron_supporting_evidence
  karyotype
  marker_feature
  misc_feature
  prediction_exon
  prediction_transcript
  protein_align_feature
  repeat_feature
  simple_feature
  transcript
);

sub usage {
	print <<USAGE_END;
USAGE:

  $0 --dbhost=ens-staging1 [--dbport=3306] \\
  \t--dbuser=ensadmin --dbpass=XXX \\
  \t--dbpattern=core

  $0 --help

  --dbpattern   Specifies a regular expression for (possibly) matching
                multiple databases.

  --help        Displays this help text.

This script will dump the current meta_coord table to a backup file in
the current directory.  Then it will update the meta_coord table for the
data in the following tables:

USAGE_END

	print( "\t", join( "\n\t", @table_names ), "\n" );

}

if ( scalar(@ARGV) == 0 ) {
	usage();
	exit 0;
}

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $optsd = $cli_helper->get_dba_opts();
# process the command line with the supplied options plus a help subroutine
my $opts = $cli_helper->process_args( $optsd, \&usage );

# use the command line options to get an array of database details
# only process each database name once (to avoid duplication for multispecies dbs)
for my $db_args ( @{ $cli_helper->get_dba_args_for_opts( $opts, 0 ) } ) {

	my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(%$db_args);

	my $file =
	  $dba->dbc()->dbname() . "_" . $dba->species_id() . ".meta_coord.backup";
	my $sys_call = sprintf( "mysql "
							  . "--host=%s "
							  . "--port=%d "
							  . "--user=%s "
							  . "--pass='%s' "
							  . "--database=%s "
							  . "--skip-column-names "
							  . " --execute='SELECT * FROM meta_coord'"
							  . " > $file",
							$dba->dbc->host(),     $dba->dbc->port(),
							$dba->dbc->username(), $dba->dbc->password(),
							$dba->dbc->dbname() );
	unless ( system($sys_call) == 0 ) {
		warn "Can't dump the original meta_coord for back up "
		  . "(skipping this species)\n";
		next;
	} else {
		print "Original meta_coord table backed up in $file\n";
	}

	foreach my $table_name (@table_names) {
		print("Updating $table_name table entries... ");

		$dba->dbc()->sql_helper()->execute_update(
			-SQL =>
"DELETE mc.* FROM meta_coord mc, coord_system cs WHERE cs.coord_system_id=mc.coord_system_id AND table_name = ? AND cs.species_id=?",
			-PARAMS => [ $table_name, $dba->species_id() ] );

                my $sql = "INSERT INTO meta_coord "
			  . "SELECT '$table_name', s.coord_system_id, "
			  . "MAX( t.seq_region_end - t.seq_region_start + 1 ) "
			  . "FROM $table_name t, seq_region s, coord_system c "
			  . "WHERE t.seq_region_id = s.seq_region_id AND c.coord_system_id=s.coord_system_id AND c.species_id=?"
			  . "GROUP BY s.coord_system_id";

		$dba->dbc()->sql_helper()->execute_update(
			-SQL => $sql,
			-PARAMS => [ $dba->species_id() ] );

		print("done\n");
        }
     

	print(   "==> Done with "
		   . $dba->dbc->dbname() . "/"
		   . $dba->species_id()
		   . "\n" );
} ## end for my $db_args ( @{ $cli_helper...})

print("==> All done.\n");

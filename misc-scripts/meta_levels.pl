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


# Populate meta table with (e.g.) genebuild.level = toplevel if all genes are
# top level. Using v41 API code this can speed fetching & dumping greatly.
#

use strict;
use warnings;

use DBI;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::CliHelper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

our @feature_types =
  qw[gene transcript exon repeat_feature dna_align_feature protein_align_feature simple_feature prediction_transcript prediction_exon];

# get the basic options for connecting to a database server
my $optsd = $cli_helper->get_dba_opts();
# add the print option
push(@{$optsd},"print");
# process the command line with the supplied options plus a help subroutine
my $opts = $cli_helper->process_args($optsd,\&usage);

# use the command line options to get an array of database details
for my $db_args (@{$cli_helper->get_dba_args_for_opts($opts)}) {
    # use the args to create a DBA
    my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$db_args});
    if($dba) {
    	process_dba($dba,$opts->{print});
    }
}

sub process_dba {

	my ($dba,$print) = @_;
	my $ma = $dba->get_MetaContainer();

	my @inserted;
	my @not_inserted;

	foreach my $type (@feature_types) {

		delete_existing( $ma, $type ) if ( !$print );

		if ( can_use_key( $dba, $type ) ) {

			insert_key( $ma, $type ) if ( !$print );
			push @inserted, $type;

		} else {

			push @not_inserted, $type;

		}

	}

	print $dba->dbc()->dbname()
	  . " (species_id "
	  . $dba->species_id()
	  . ") inserted keys for "
	  . join( ", ", @inserted ) . ".\n"
	  if (@inserted);
	print "Did not insert keys for " . join( ", ", @not_inserted ) . ".\n"
	  if (@not_inserted);
} ## end sub process_dba

#------------------------------------------------------------------------------

sub delete_existing {

	my ( $ma, $type ) = @_;

	$ma->delete_key( $type . "build.level" );

}

#------------------------------------------------------------------------------

sub can_use_key {

	my ( $dba, $type ) = @_;

 # compare total count of typewith number of toplevel type, if they're the same,
 # then we can use the key

	my $sth = $dba->dbc()->prepare("SELECT COUNT(*) FROM $type");
	$sth->execute();
	my $total = ( $sth->fetchrow_array() )[0];

	$sth =
	  $dba->dbc()
	  ->prepare(   "SELECT COUNT(*) "
				 . "FROM $type t, seq_region_attrib sra, attrib_type at "
				 . "WHERE t.seq_region_id=sra.seq_region_id "
				 . "AND sra.attrib_type_id=at.attrib_type_id "
				 . "AND at.code='toplevel'" );
	$sth->execute();
	my $toplevel = ( $sth->fetchrow_array() )[0];

	if ( $toplevel > 0 ) {
		return $total == $toplevel;
	}
} ## end sub can_use_key

#------------------------------------------------------------------------------

sub insert_key {
	my ( $ma, $type ) = @_;
	$ma->store_key_value( $type . "build.level", "toplevel" );
}

#------------------------------------------------------------------------------

sub usage {
	print <<EOF; exit(0);

Populate meta table with (e.g.) genebuild.level = toplevel if all genes
are top level. Using v41 API code this can speed fetching and dumping
greatly.

Usage: perl $0 <options>

  --host|dbhost     Database host to connect to.

  --port|dbport     Database port to connect to (default is 3306).

  --dbpattern       Database name regexp

  --user|dbuser     Database username.

  --pass|dbpass     Password for user.

  --print           Just print, don't insert or delete keys.

  --help            This message.

  --verbose         Prints out the name of the database
                    which is going to be patched.



   If you like to patch just a single database you can also use this
   command line :

   perl meta_levels.pl --host ... --user ... --pass ... --port ... --dbname ...
EOF
} ## end sub usage

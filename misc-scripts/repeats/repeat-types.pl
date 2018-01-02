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

# Repeat classification script
#
# This script is used to do the repeat classification for web display
# on newer v32 databases.
#

use strict;
use warnings;

use Bio::EnsEMBL::Utils::CliHelper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $optsd = $cli_helper->get_dba_opts();
# add the print option
push( @{$optsd}, "print|p" );
# process the command line with the supplied options plus a help subroutine
my $opts = $cli_helper->process_args( $optsd, \&usage );

# use the command line options to get an array of database details
for my $db_args ( @{ $cli_helper->get_dba_args_for_opts($opts) } ) {

	# use the args to create a DBA
	my $dba    = Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{$db_args} );
	my $helper = $dba->dbc()->sql_helper();
	print STDOUT "Processing species "
	  . $dba->species_id()
	  . " from database "
	  . $dba->dbc()->dbname()
	  . " on server "
	  . $dba->dbc()->host() . "\n";
	print STDERR "  Setting repeat types\n";

	my %mappings = ( 'Low_Comp%'  => 'Low complexity regions',
					 'LINE%'      => 'Type I Transposons/LINE',
					 'SINE%'      => 'Type I Transposons/SINE',
					 'DNA%'       => 'Type II Transposons',
					 'LTR%'       => 'LTRs',
					 'Other%'     => 'Other repeats',
					 'Satelli%'   => 'Satellite repeats',
					 'Simple%'    => 'Simple repeats',
					 'Other%'     => 'Other repeats',
					 'Tandem%'    => 'Tandem repeats',
					 'TRF%'       => 'Tandem repeats',
					 'Waterman'   => 'Waterman',
					 'Recon'      => 'Recon',
					 'Tet_repeat' => 'Tetraodon repeats',
					 'MaskRegion' => 'Mask region',
					 'dust%'      => 'Dust',
					 'Unknown%'   => 'Unknown',
					 '%RNA'       => 'RNA repeats', );
	foreach ( keys %mappings ) {
		$helper->execute_update( -SQL =>
qq(update repeat_consensus set repeat_type = '$mappings{$_}' where repeat_class like '$_')
		);
	}

	# type all remaining repeats as unknown
	$helper->execute_update( -SQL =>
qq(update repeat_consensus set repeat_type = 'Unknown' where repeat_type = '')
	);
	$helper->execute_update( -SQL =>
qq(update repeat_consensus set repeat_type = 'Unknown' where repeat_type is null)
	);
} ## end for my $db_args ( @{ $cli_helper...})

print STDERR "All done.\n";

sub usage {
	print STDERR <<EOF

This program classifies the repeats stored in a core database into some
somewhat sensible categories.  It does this through a combination of a
repeat.txt file extracted from RepeatMasker repeat libraries and through
some simple pattern matching of the repeat names.

usage: perl repeat-types.pl  [-user <user>] [-port <port>] [-pass <pass>]
               -host <host> -dbpattern <regexp>

example: perl repeat-types.pl -user ensadmin -pass secret -host ecs1g \\
             -port 3306 -dbpattern '^homo_sapiens_(core|vega)_20_34c$'

EOF
	  ;
	exit;
}

#!/usr/bin/env perl

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Getopt::Long;
use Bio::EnsEMBL::Utils::CliHelper;

my @table_names = qw(
  assembly_exception
  gene
  exon
  density_feature
  ditag_feature
  dna_align_feature
  karyotype
  marker_feature
  misc_feature
  qtl_feature
  prediction_exon
  prediction_transcript
  protein_align_feature
  repeat_feature
  simple_feature
  splicing_event
  transcript
);

sub usage {
	my ($fail) = @_;
	my $usage = qq(
  $0 --host ens-staging --port 3306 --user ensadmin \\
    --pass XXX --dbpattern core

  [--help] displays this menu.

This script will dump the current meta_coord table in the latest
homo_sapiens_core.meta_coord file.  Then it will update the meta_coord
table for all the following table names one by one

  assembly_exception
  gene
  exon
  density_feature
  ditag
  ditag_feature
  dna_align_feature
  karyotype
  marker_feature
  misc_feature
  qtl_feature
  prediction_exon
  prediction_transcript
  protein_align_feature
  repeat_feature
  simple_feature
  splicing_event
  transcript
  );
	exit $fail;
} ## end sub usage

use Bio::EnsEMBL::Utils::CliHelper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $optsd = $cli_helper->get_dba_opts();
# process the command line with the supplied options plus a help subroutine
my $opts = $cli_helper->process_args( $optsd, \&usage );

# use the command line options to get an array of database details
# only process each database name once (to avoid duplication for multispecies dbs)
for my $db_args ( @{ $cli_helper->get_dba_args_for_opts( $opts, 1 ) } ) {

	my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(%$db_args);
	my $file =
	  $dba->dbc()->dbname() . "_" . $dba->species_id() . ".meta_coord.backup";
	my $sys_call = sprintf( "mysql -h%s -P%d -u%s -p'%s' %s > $file",
							$dba->dbc->host(), $dba->dbc->port(),
							$dba->dbc->user(), $dba->dbc->pass(),
							$dba->dbc->dbname() );

	unless ( system($sys_call) == 0 ) {
		print STDERR "Can't dump the original meta_coord for back up\n";
		exit 1;
	} else {
		print STDERR "Original meta_coord table backed up in $file\n";
	}

	foreach my $table_name (@table_names) {

		print STDERR "Updating $table_name table entries...";
		$dba->dbc()->helper()->execute_update(
			-SQL =>
"DELETE mc.* FROM meta_coord mc, coord_system cs WHERE cs.coord_system_id=mc.coord_system_id AND table_name = ? AND cs.species_id=?"
			,
			-PARAMS => [ $table_name, $dba->species_id() ] );

		$dba->dbc()->helper()->execute_update(
			-SQL => "INSERT INTO meta_coord "
			  . "SELECT '$table_name', s.coord_system_id, "
			  . "MAX( t.seq_region_end - t.seq_region_start + 1 ) "
			  . "FROM $table_name t, seq_region s, coord_system c "
			  . "WHERE t.seq_region_id = s.seq_region_id AND c.coord_system_id=s.coord_system_id AND c.species_id=?"
			  . "GROUP BY s.coord_system_id",
			-PARAMS => [ $dba->species_id() ] );
		print STDERR "Done\n";
	}
} ## end for my $db_args ( @{ $cli_helper...})

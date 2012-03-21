#!/usr/bin/env perl
# Calculate per-gene GC content and store as gene attributes

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
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
	my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{$db_args} );

	print STDOUT "Processing species "
	  . $dba->species_id()
	  . " from database "
	  . $dba->dbc()->dbname()
	  . " on server "
	  . $dba->dbc()->host() . "\n";
	delete_existing($dba) unless ( $opts->{print} );

	print STDOUT "Calculating Gene GC attributes\n";
	my $attribute_adaptor = $dba->get_AttributeAdaptor();
	my $total_count       = 0;
	for my $gene ( @{ $dba->get_GeneAdaptor()->fetch_all() } ) {
		my $gc = $gene->feature_Slice()->get_base_count->{'%gc'};
		if ( !$opts->{print} ) {
	 # attribute types need to be propagated from production database to all dbs
	 # if the type exists it won't be overwritten
			my $attribute =
			  Bio::EnsEMBL::Attribute->new(
						  -CODE        => 'GeneGC',
						  -NAME        => 'Gene GC',
						  -DESCRIPTION => 'Percentage GC content for this gene',
						  -VALUE       => $gc );
			my @attributes = ($attribute);
			$attribute_adaptor->store_on_Gene( $gene->dbID, \@attributes );
			$total_count++;
		} else {
			print $gene->stable_id() . " " . $gc . "\n";
		}
	}

	if ( !$opts->{print} ) {
		print STDOUT "Written $total_count 'GeneGC' gene attributes to species "
		  . $dba->species_id()
		  . " from database "
		  . $dba->dbc()->dbname()
		  . " on server "
		  . $dba->dbc()->host() . "\n";
	}

} ## end for my $db_args ( @{ $cli_helper...})

# ----------------------------------------------------------------------

sub delete_existing {

	my $dba = shift;

	print STDOUT "Deleting existing 'GeneGC' gene attributes\n";
	$dba->dbc()->sql_helper()->execute_update(
		-SQL =>
q/DELETE ga FROM gene_attrib ga, attrib_type at, gene g, seq_region s, coord_system c 
WHERE at.attrib_type_id=ga.attrib_type_id AND at.code='GeneGC'
AND ga.gene_id=g.gene_id AND g.seq_region_id=s.seq_region_id 
AND c.coord_system_id=s.coord_system_id AND c.species_id=?
/,
		PARAMS => [ $dba->species_id() ] );

	return;
}

sub usage {
	my $indent = ' ' x length($0);
	print <<EOF; exit(0);

What does it do?

This script calculates per-gene GC content and stores it as gene attributes.
It deletes existing Gene GC attributes. Then fetches all genes in the
core db and calculates the %gc for each gene.

Input data: dna sequence 
Output tables: gene_attrib 


When to run it in the release cycle?

It can be run whenever the genes and sequence are stable, i.e. any time after 
the genebuild handover, but before the handover to Mart.


Which databases to run it on?

It needs to be run across all core databases for every release.


How long does it take?

It takes a total of about 10 hours to run for all core databases in normal queue,


Usage:

  $0 -h host [-port port] -u user -p password \\
  $indent -pattern pattern [-print] \\
  $indent [-help]  \\

  -h|host              Database host to connect to

  -port                Database port to connect to (default 3306)

  -u|user              Database username

  -p|pass              Password for user

  -pattern             Database name regexp

  -print               Just print, don't insert or delete attributes

  -help                This message


EOF

} ## end sub usage

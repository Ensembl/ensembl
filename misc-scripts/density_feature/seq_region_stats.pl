#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::CliHelper;
use Data::Dumper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();
my $optsd =
  [ @{ $cli_helper->get_dba_opts() }, @{ $cli_helper->get_dba_opts('m') } ];
# add the print option
push( @{$optsd}, "stats|s:s" );
my $opts = $cli_helper->process_args( $optsd, \&usage );

usage() if ( !$opts->{stats} || $opts->{stats} !~ /^(gene|snp)$/ );

if ( !defined $opts->{mhost} ) {
	( $opts->{mhost}, $opts->{mport}, $opts->{mdbname}, $opts->{muser} ) =
	  ( 'ens-staging1', '3306', 'ensembl_production', 'ensro' );
}

my ($prod_dba) = @{ $cli_helper->get_dbas_for_opts( $opts, 1, 'm' ) };

if ( !defined $prod_dba ) {
	usage();
}

my %attrib_codes = %{
	$prod_dba->dbc()->sql_helper()->execute_into_hash(
		-SQL =>
"select distinct b.name, code from biotype b join attrib_type using(attrib_type_id) where is_current = 1 and db_type like '%core%' and object_type = 'gene' order by b.name"
	),
	,
	{ Columns => [ 1, 2 ] } };

#add known and novel protein coding attrib types
$attrib_codes{'known protein_coding'} = 'GeneNo_knwCod';
$attrib_codes{'novel protein_coding'} = 'GeneNo_novCod';

my $genestats = 1 if ( $opts->{stats} eq 'gene' );
my $snpstats  = 1 if ( $opts->{stats} eq 'snp' );
for my $db_args ( @{ $cli_helper->get_dba_args_for_opts($opts) } ) {

	my $db     = new Bio::EnsEMBL::DBSQL::DBAdaptor(%$db_args);
	my $dbname = $db->dbc()->dbname();

	my $total_count = 0;
	# delete old attributes before starting
	if ($genestats) {
		foreach my $code ( values %attrib_codes ) {
			$db->dbc()->sql_helper()->execute_update(
				-SQL =>
"DELETE sa FROM seq_region_attrib sa, attrib_type at, seq_region s, coord_system c WHERE s.seq_region_id=sa.seq_region_id AND c.coord_system_id=s.coord_system_id AND at.attrib_type_id=sa.attrib_type_id AND at.code=? AND c.species_id=?",
				-PARAMS => [ $code, $db->species_id() ] );
		}
	}

	if ($snpstats) {
		$db->dbc()->sql_helper()->execute_update(
			-SQL =>
"DELETE sa FROM seq_region_attrib sa, attrib_type at, seq_region s, coord_system c WHERE s.seq_region_id=sa.seq_region_id AND c.coord_system_id=s.coord_system_id AND at.attrib_type_id=sa.attrib_type_id AND at.code=? AND c.species_id=?",
			-PARAMS => [ "SNPCount", $db->species_id() ] );
	}

	#
	# Only run on database with genes
	#

	my $genes_present;

	if ($genestats) {
		my $gene_count =
		  $db->dbc()->sql_helper()->execute_single_result(
			-SQL =>
"select count(*) from gene join seq_region using (seq_region_id) join coord_system using (coord_system_id) where species_id=?",
			-PARAMs => [ $db->species_id() ] );
		$genes_present = ($gene_count) ? 1 : 0;
	} else {
		$genes_present = 0;
	}

	#
	# and seq_regions
	#
	my $seq_region_count =
	  $db->dbc()->sql_helper()->execute_single_result(
		-SQL =>
"select count(*) from seq_region join coord_system using (coord_system_id) where species_id=?",
		-PARAMS => [ $db->species_id() ] );
	if ( !$seq_region_count ) {
		print STDERR "No seq_regions for $dbname\n";
		exit();
	}

	my $snps_present;
	my $snp_db;

	if ($snpstats) {
		$snp_db = variation_attach($db);
		if ( defined $snp_db ) { $snps_present = 1; }
	}

	my $slice_adaptor  = $db->get_SliceAdaptor();
	my $attrib_adaptor = $db->get_AttributeAdaptor();

	# Do not include non-reference sequences ie. haplotypes for human
	#my $top_slices = $slice_adaptor->fetch_all( "toplevel" , undef, 1);
	my $top_slices = $slice_adaptor->fetch_all("toplevel");

	while ( my $slice = shift( @{$top_slices} ) ) {
	#    print STDERR "Processing seq_region ", $slice->seq_region_name(), "\n";

		my @attribs;

		if ($genes_present) {

			my %attrib_counts;
			my %counts;

			my $genes = $slice->get_all_Genes();

			while ( my $gene = shift( @{$genes} ) ) {

				my $biotype = $gene->biotype();
				if ( $biotype =~ /coding/i && $biotype !~ /non_/i ) {
					if ( $gene->is_known() ) {
						$biotype = "known " . $biotype;
					} else {
						$biotype = "novel " . $biotype;
					}
				}

				$counts{$biotype}++;

			}

			for my $biotype ( keys %counts ) {
				my $attrib_code = $attrib_codes{$biotype};
				if ( !$attrib_code ) {
					print STDERR
					  "Unspecified biotype \"$biotype\" in database $dbname.\n";
					next;
				}

				$attrib_counts{$attrib_code} += $counts{$biotype};

			}

			foreach my $attrib_code ( keys %attrib_counts ) {
				my $attrib_code_desc = $attrib_code;
				$attrib_code_desc =~ s/GeneNo_//;
				push @attribs,
				  Bio::EnsEMBL::Attribute->new(
					 -NAME        => $attrib_code_desc . ' Gene Count',
					 -CODE        => $attrib_code,
					 -VALUE       => $attrib_counts{$attrib_code},
					 -DESCRIPTION => 'Number of ' . $attrib_code_desc . ' Genes'
				  );

			}

		} ## end if ($genes_present)

		if ($snps_present) {
			my $count =
			  $snp_db->dbc()->sql_helper()->execute_single_result(
				-SQL =>
"SELECT COUNT(*) FROM variation_feature WHERE seq_region_id = ?",
				-PARAMS => [ $slice->get_seq_region_id ] );

			push @attribs,
			  Bio::EnsEMBL::Attribute->new(
										  -NAME        => 'SNP Count',
										  -CODE        => 'SNPCount',
										  -VALUE       => $count,
										  -DESCRIPTION => 'Total Number of SNPs'
			  );
		}

		$attrib_adaptor->store_on_Slice( $slice, \@attribs );
		my $slice_attrib_count = @attribs;
		$total_count += $slice_attrib_count;
		#  print_chromo_stats([$slice]);
	} ## end while ( my $slice = shift...)

	print STDOUT
"Written $total_count seq reqion attributes to database $dbname on server "
	  . $db->dbc()->host() . "\n";

} ## end for my $db_args ( @{ $cli_helper...})

sub print_chromo_stats {
	my $chromosomes = shift;

	foreach my $chr (@$chromosomes) {
		print "\nchromosome: ", $chr->seq_region_name(), "\n";
		foreach my $attrib ( @{ $chr->get_all_Attributes() } ) {
			print "  ", $attrib->name(), ": ", $attrib->value(), "\n";
		}
	}
}

#
# tries to attach variation database.
#

sub variation_attach {
	my $db = shift;

	my $core_db_name;
	$core_db_name = $db->dbc->dbname();
	if ( $core_db_name !~ /_core_/ ) {
		return 0;
	}
	#
	# get a lost of all databases on that server
	#
	my $sth = $db->dbc->prepare("show databases");
	$sth->execute();
	my $all_db_names = $sth->fetchall_arrayref();
	my %all_db_names = map { ( $_->[0], 1 ) } @$all_db_names;
	my $snp_db_name  = $core_db_name;
	$snp_db_name =~ s/_core_/_variation_/;

	if ( !exists $all_db_names{$snp_db_name} ) {
		return 0;
	}

	# this should register the dbadaptor with the Registry
	my $snp_db =
	  Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
												-host => $db->dbc()->host(),
												-user => $db->dbc()->username(),
												-pass => $db->dbc()->password(),
												-port => $db->dbc()->port(),
												-dbname  => $snp_db_name,
												-group   => "variation",
												-species => "DEFAULT" );

	return $snp_db;
} ## end sub variation_attach

sub usage {
	my $indent = ' ' x length($0);
	print <<EOF; exit(0);


For each toplevel slice, count the number of genes for each biotype
(gene stats) or count the number of SNPs (snp stats).

gene stats
What does it do?

Deletes all seq_region_attrib that have attrib_type code 
with a prefix 'GeneNo_'. All toplevel slices are fetched.

Input data: dna seqence, genes, xrefs, xref projections 
Output tables: seq_region_attrib (attrib_type code with prefix 'GeneNo')


When to run it in the release cycle?

After core have finished xref projections


Which databases to run it on?

Run on all core databases (including otherfeatures, cdna etc) for each release.


How long does it take?

It takes about 10 mins to run for a database in normal queue,


snp stats

What does it do?

Deletes out all seq_region_attrib that have attrib_type code of 'SNPCount'. 
Attach variation db if exists. All toplevel slices are fetched. 
For each slice, count the number of SNPs.

This option requires ensembl-variation in perl5lib.

Input data: top level seq regions, variation db
Output tables: seq_region_attrib (attrib_type code with prefix 'SNPCount')


When to run it in the release cycle?

When variation dbs have been handed over


Which databases to run it on?

Run on core databases only for new species or if the assembly changed, 
or if the variation positions have changed in the corresponding variation db.


How long does it take?

It takes about 20 mins to run for a database in normal queue.



Usage: 

  $0 -h host [-port port] -u user -p password \\
  $indent -d database | -pattern pattern \\
  $indent [-mhost ensembl_production host] [-mport ensembl_production host] [-muser ensembl_production host] \\
  $indent -s gene | snp  \\
  $indent [-help]  \\

  -h|host             Database host to connect to

  -port               Database port to connect to (default 3306)

  -u|user             Database username

  -p|pass             Password for user

  -d|dbname           Database name

  -pattern            Database name regexp

  -mhost              ensembl_production database host to connect to

  -mport              ensembl_production database port to connect to

  -muser              ensembl_production database username

  -s|stats            'gene' or 'snp'

  -help               This message


EOF

} ## end sub usage

#!/usr/bin/env perl

#
# Calculate the GC content for top level seq_regions
#   small regions 500bp to be able to display on contigview
#   big regions genomesize / 4000 for 4000 features on the genome

use strict;
use warnings;
#use dbi;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DensityTypeAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Utils::CliHelper;
use Data::Dumper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();
my $optsd =
  [ @{ $cli_helper->get_dba_opts() }, @{ $cli_helper->get_dba_opts('m') } ];
push( @{$optsd}, "print" );
my $opts = $cli_helper->process_args( $optsd, \&usage );
if ( !defined $opts->{mdbname} && !defined $opts->{mdbpattern} ) {
	$opts->{mdbname} = 'ensembl_production';
}

if ( !defined $opts->{mhost} ) {
	( $opts->{mhost}, $opts->{mport}, $opts->{mdbname}, $opts->{muser} ) =
	  ( 'ens-staging1', '3306', 'ensembl_production', 'ensro' );
}

my ($prod_dba) = @{ $cli_helper->get_dbas_for_opts( $opts, 1, 'm' ) };

if ( !defined $prod_dba ) {
	usage();
}

for my $db_args ( @{ $cli_helper->get_dba_args_for_opts($opts) } ) {
	print STDOUT "Connecting to " . $db_args->{-DBNAME} . "\n";

	my $db     = new Bio::EnsEMBL::DBSQL::DBAdaptor(%$db_args);
	my $dbname = $db->dbc()->dbname();

	my $helper = $db->dbc()->sql_helper();

	my $bin_count  = 150;
	my $max_slices = 100;

	#
	# Check wether the script should run on given database
	#
	my $dna_count = $helper->execute_single_result(
		-SQL => q/select count(*) from dna join seq_region 
using (seq_region_id) join coord_system using (coord_system_id) where species_id=?/,
		-PARAMS => [ $db->species_id() ] );

	if ( !$dna_count ) {
		print STDERR "No dna, no gc content for $dbname.\n";
		exit();
	}
	print STDOUT "Deleting old PercentGC features\n";
	$helper->execute_update(
		-SQL => qq(
DELETE df, dt
FROM density_feature df, density_type dt, analysis a,
seq_region sr, coord_system c
WHERE 
sr.seq_region_id=df.seq_region_id 
AND sr.coord_system_id=c.coord_system_id
AND c.species_id=?
AND a.analysis_id=dt.analysis_id
AND dt.density_type_id=df.density_type_id
AND a.logic_name='percentgc'),
		-PARAMS => [ $db->species_id() ] );

	#
	# Get the adaptors needed;
	#

	my $slice_adaptor = $db->get_SliceAdaptor();
	my $dfa           = $db->get_DensityFeatureAdaptor();
	my $dta           = $db->get_DensityTypeAdaptor();
	my $aa            = $db->get_AnalysisAdaptor();

	# Sort slices by coordinate system rank, then by length.
	my @sorted_slices =
	  sort( {      $a->coord_system()->rank() <=> $b->coord_system()->rank()
				|| $b->seq_region_length() <=> $a->seq_region_length()
	  } @{ $slice_adaptor->fetch_all('toplevel') } );

	my $analysis = $aa->fetch_by_logic_name('percentgc');

	if ( !defined($analysis) ) {

		my ( $display_label, $description ) =
		  @{$prod_dba->dbc()->sql_helper()->execute(
				-SQL =>
"select distinct display_label, description from analysis_description where is_current = 1 and logic_name = 'percentgc'"
			) };

		$analysis =
		  new Bio::EnsEMBL::Analysis( -program       => "percent_gc_calc.pl",
									  -database      => "ensembl",
									  -gff_source    => "percent_gc_calc.pl",
									  -gff_feature   => "density",
									  -logic_name    => "percentgc",
									  -description   => $description,
									  -display_label => $display_label,
									  -displayable   => 1 );

		$aa->store($analysis);

	} else {

		my $support = 'Bio::EnsEMBL::Utils::ConversionSupport';
		$analysis->created( $support->date() );
		$aa->update($analysis);

	}

	#
	# Create new density type.
	#

	my $density_type =
	  Bio::EnsEMBL::DensityType->new( -analysis        => $analysis,
									  -region_features => $bin_count,
									  -value_type      => 'ratio' );

	$dta->store($density_type);

	my ( $current_start, $current_end, $current );
	my $slice_count = 0;
	my $block_size;
	my $total_count;
	while ( my $slice = shift @sorted_slices ) {

		$block_size = $slice->length()/$bin_count;

		my @density_features = ();

		print STDOUT "GC percentage for "
		  . $slice->seq_region_name()
		  . " with block size $block_size\n";

		$current_end = 0;
		$current     = 0;

		while ( $current_end < $slice->length ) {

			$current += $block_size;
			$current_start = $current_end + 1;
			$current_end   = int( $current + 1 );

			if ( $current_end < $current_start ) {
				$current_end = $current_start;
			}

			if ( $current_end > $slice->end() ) {
				$current_end = $slice->end();
			}

			my $sub_slice = $slice->sub_Slice( $current_start, $current_end );

			my $gc = $sub_slice->get_base_count()->{'%gc'};
			my $df =
			  Bio::EnsEMBL::DensityFeature->new( -seq_region => $slice,
												 -start      => $current_start,
												 -end        => $current_end,
												 -density_type => $density_type,
												 -density_value => $gc );
#print join ("\t", $slice, $current_start, $current_end, $density_type, $gc, "\n");
			push( @density_features, $df );
			if ( $gc > 0 ) {
                                $dfa->store($df);
				#density features with value = 0 are not stored
				$total_count++;
			}
		} ## end while ( $current_end < $slice...)

		last if ( $slice_count++ > $max_slices );
	} ## end while ( my $slice = shift...)

	print STDOUT "Created $total_count percent gc density features\n";
} ## end for my $db_args ( @{ $cli_helper...})

#
# helper to draw an ascii representation of the density features
#
sub print_features {
	my $features = shift;

	return if ( !@$features );

	my $sum    = 0;
	my $length = 0;
	#  my $type = $features->[0]->{'density_type'}->value_type();

	print("\n");
	my $max = 0;
	foreach my $f (@$features) {
		if ( $f->density_value() > $max ) {
			$max = $f->density_value();
		}
	}
	foreach my $f (@$features) {
		my $i = 1;
		for ( ; $i < ( $f->density_value()/$max )*40; $i++ ) {
			print "*";
		}
		for ( my $j = $i; $j < 40; $j++ ) {
			print " ";
		}
		print "  " . $f->density_value() . "\t" . $f->start() . "\n";
	}
	#  my $avg = undef;
	#  $avg = $sum/$length if($length < 0);
	#  print("Sum=$sum, Length=$length, Avg/Base=$sum");
} ## end sub print_features

sub usage {
	my $indent = ' ' x length($0);
	print <<EOF; exit(0);


What does it do?

First, it needs the dna table to be populated. It then deletes
all PercentGC entries from the analysis, density_type and 
density_feature tables. All toplevel slices are fetched and sorted
from longest to shortest. Each slice is divided into 150 bins
(blocks). For each of the blocks or sub-slices, the %gc is
calculated.

Input data: dna sequence, top level seq regions 
Output tables: updates analysis creation date, 
               density_type, density_feature


When to run it in the release cycle?

It can be run after genebuilders have finished their Xrefs 
(script not affected by projected Xrefs).


Which databases to run it on?

Run on core databases for new species or if one of the following changed:
  - dna sequence
  - assembly


How long does it take?

It takes about 15 mins to run for a database in normal queue,


Usage:

  $0 -h host [-port port] -u user -p password \\
  $indent -d database \\
  $indent [-mhost ensembl_production host] [-mport ensembl_production host] [-muser ensembl_production host] \\
  $indent [-help]  \\

  -h|host              Database host to connect to

  -port                Database port to connect to (default 3306)

  -u|user              Database username

  -p|pass              Password for user

  -pattern             Database name regexp

  -mhost              ensembl_production database host to connect to

  -mport              ensembl_production database port to connect to

  -muser              ensembl_production database username

  -help                This message


EOF

} ## end sub usage


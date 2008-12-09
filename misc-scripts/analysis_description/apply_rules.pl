#!/opt/local/bin/perl
###!/usr/local/ensembl/bin/perl

# POD documentation - main docs before the code

=pod

=head1 NAME

  apply_rules.pl

=head1 SYNOPSIS

 script applies additonal rules to the analysis description table posterior to
 description loading

=head1 DESCRIPTION

 The rules (for details talk to Steve T.):
   1) For each species, dna_align_features are displayed according to which 
      database they can be retrieved from:
      - if they are present in both otherfeatures and core, then only the ones 
        from otherfeatures are to be displayed, ie we need to set the displayable 
        entry to 0 in core
	  - if they are only present in one of these two databases then show them 
        from that source

   2) Logic names of human_cdna and mouse_cdna are slightly different in that 
      they should be switched off (ie set the displayable entry to 0) in human 
      and mouse respectively - they are superceded by the cDNA update features 
      in the cDNA databases

   3) These cDNA_update features should have a display label of 'Mouse cDNA' in 
      the mouse_cdna database, and 'Human cDNA' in the human_cdna database

   4) All align_features from vega databases should be switched off

   5) Genes with a logic_name of 'havana' have a display label of 'Vega gene' in 
      mouse and 'Vega Havana gene' in human.

 
=head1 OPTIONS

     Database options

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname      For RDBs, what name to connect to (dbname= in locator)
    -dbuser      For RDBs, what username to connect as (dbuser= in locator)
    -dbpass      For RDBs, what password to use (dbpass= in locator)
    -file        Path to file containing descriptions. The file 
                   analysis.descriptions in this directory can be used and is 
                   supposed to be the reference file
	-update      Perform actual updates of analyses
    -help print out documentation

=head1 EXAMPLES


=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

$! = 1;

my ($dsn, $dbh, $update);
# Analysis adaptors
my ($caa, $ofaa, $cdnaaa, $vegaaa);

my $dbhost = '';
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;
my $help = 0;
my $version = 52;
my $file = 'analysis.descriptions';

&GetOptions (
	'host|dbhost=s'       => \$dbhost,
	'port|dbport=s'       => \$dbport,
	'user|dbuser=s'       => \$dbuser,
	'pass|dbpass=s'       => \$dbpass,
	'dbname=s'            => \$dbname,
	'version=s'           => \$version,
	'file|descriptions=s' => \$file,
	'update'              => \$update,
	'h|help!'             => \$help
	);

if(!$dbhost){
  print ("Need to pass in -dbhost $dbhost and -dbname $dbname\n");
  $help = 1;
}

if($help){
  usage();
}

my %reference;
if ($file) {

	open(FH, $file) 
		or throw("Failed to open reference file '$file': $@");

	while (<FH>) {
		
		chomp;
		next if m/^\#/;   # skip comments
		next if m/^$/;    # and blank lines
		next if m/^\s+$/; # and whitespace-only lines
    
		my ($nr, $logic_name, $description, $display_label, $displayable, $web_data) = split(/\t/);
		#print join("\t", $logic_name, $description, $display_label, $displayable, $web_data), "\n";

		warn ("Displayable flag for analysis '$logic_name' has to be either 0 or 1, but not '$displayable'!")
			unless ($displayable =~ m/^[01]$/);
		
		$reference{lc($logic_name)} = {
			nr            => "$nr",
			description   => "$description",
			display_lable => "$display_label", 
			displayable   => "$displayable", 
			web_data      => "$web_data"
		};
		
	}

	close FH;
	
} else {

	throw("Need to pass reference file with analysis descriptions!");

}

$dsn = "DBI:mysql:host=" . $dbhost . ";port=" . $dbport;
eval{ 
	$dbh = DBI->connect($dsn, $dbuser, $dbpass, 
						{'RaiseError' => 1,
						 'PrintError' => 0});
};


# get core database(s);
my $pat = defined $dbname ? $dbname : "%core_$version%";
my $sql = "show databases like '$pat'";
my $cdbs  = $dbh->selectcol_arrayref($sql);

foreach my $cdb (@$cdbs) {

	(my $species = $cdb) =~ s/(.+)_core_${version}_\d+[a-z]$/$1/;
	#print Dumper $species;

	my $cdba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
		-host    => $dbhost,
		-user    => $dbuser,
		-dbname  => $cdb,
		-pass    => $dbpass,
		-port    => $dbport,
		-species => $species
		);
	#print Dumper $cdba;

	$caa = $cdba->get_AnalysisAdaptor();

	### implements rule 1. ###
	
	my $cdaf_logic_names = get_af_logic_names($cdba, 'dna');
	#print Dumper $cdaf_logic_names;

	(my $ofdb = $cdb) =~ s/_core_/_otherfeatures_/;
	$sql = "show databases like '$ofdb'";
	my $ofdbs = $dbh->selectcol_arrayref($sql);

	if (scalar(@$ofdbs) == 0) {
		
		print ("No otherfeatures db for " . $cdb . "! Setting all displayable entires to 1\n");
		my $daf_logic_names = get_af_logic_names($cdba, 'dna');

		map { update_analysis($caa, $_, 1) } @$daf_logic_names;

	} else {

		print ("Both core and otherfeatures dbs exist. Need to analyse dna_align_features ...\n");

		my $ofdba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
			-host    => $dbhost,
			-user    => $dbuser,
			-dbname  => $ofdbs->[0],
			-pass    => $dbpass,
			-port    => $dbport,
			-species => $species,
			-group   => 'otherfeatures'
			);
		#print Dumper $ofdba;

		$ofaa = $ofdba->get_AnalysisAdaptor();
		
		my $ofdaf_logic_names = get_af_logic_names($ofdba, 'dna');
		#print Dumper $ofdaf_logic_names;

		my %daf_logic_names;
		map {$daf_logic_names{lc($_)}++} (@$cdaf_logic_names, @$ofdaf_logic_names);

		foreach my $ln (@$ofdaf_logic_names) {

			if ($daf_logic_names{lc($ln)} == 2) {

				print("<$ln> exists in both, setting displayable 0 for core and 1 for otherfeatures\n");
				update_analysis($caa, $ln, 0);
				update_analysis($ofaa, $ln, 1);

			} else {
				
				print("<$ln> exists only in otherfeatures, setting displayable according to reference file\n");
				update_analysis($ofaa, $ln, $reference{lc($ln)}{displayable});
	
			}
			delete $daf_logic_names{lc($ln)};

		}

		foreach my $ln (keys %daf_logic_names) {

				print("<$ln> exists only in core, setting displayable according to reference file\n");
				update_analysis($caa, $ln, $reference{lc($ln)}{displayable});

		}

	}
	

	if ($species =~ m/^(homo_sapiens|mus_musculus)/) {


		### implements rule 5. ###

		my $display_label = ($species eq 'homo_sapiens') ? 'Vega Havana gene' : 'Vega gene';
		print "<$cdb> Updating display_label for logic_name havana to '$display_label'\n";
		update_analysis($caa, 'havana', undef, $display_label);

		### implements rule 2. ans 3. ###

		(my $cdnadb = $cdb) =~ s/_core_/_cdna_/;
		my $cdnadba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
			-host    => $dbhost,
			-user    => $dbuser,
			-dbname  => $cdnadb,
			-pass    => $dbpass,
			-port    => $dbport,
			-species => $species,
			-group   => 'cdna'
			);
		#print Dumper $cdnadba;
		$cdnaaa = $cdnadba->get_AnalysisAdaptor();

		my %alias = (
			'homo_sapiens' => 'Human',
			'mus_musculus' => 'Mouse'
		);

		my $ln = lc($alias{$species}).'_cdna';
		print "<$cdb> Switching off displayable for $ln\n";
		update_analysis($caa, $ln, 0);

		my $dl = $alias{$species}.' cDNA';
		print "<$cdnadb> Updating display_label for cDNA_update to '$dl'\n";
		update_analysis($cdnaaa, 'cDNA_update', 1, $dl);
			

		### implements rule 4. ###
		
		(my $vegadb = $cdb) =~ s/_core_/_vega_/;
		my $vegadba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
			-host    => $dbhost,
			-user    => $dbuser,
			-dbname  => $vegadb,
			-pass    => $dbpass,
			-port    => $dbport,
			-species => $species,
			-group   => 'vega'
			);
		#print Dumper $vegadba;
		$vegaaa = $vegadba->get_AnalysisAdaptor();
	
		my $vega_daf_logic_names = get_af_logic_names($vegadba, 'dna');
		my $vega_paf_logic_names = get_af_logic_names($vegadba, 'protein');
		
		foreach my $ln (@$vega_daf_logic_names, @$vega_paf_logic_names) {

			print "<$vegadb> Switching align_feature '$ln' displayable off\n";
			update_analysis($vegaaa, $ln, 0);
			

		}

		# 6) There are two more rules for the logic_name 'ensembl' in core databases 
		#    depending on the species:
		#	 - all but mouse, human and anopheles, use definition file.
		#	 - anopheles needs a display_label of 'VectorBase gene'
		#	 - mouse and human need a different web_data column:
		#	   '{'caption' => 'Ensembl/Havana gene','name' => 'Merged Ensembl and Havana Genes',
		#        'label_key' => '[text_label] [display_label]',
		#        'default' => {'contigviewbottom' => 'transcript_label',
		#        'contigviewtop' => 'gene_label',
		#        'cytoview' => 'gene_label'},'key' => 'ensembl'}'

		my $web_data = "'{'caption' => 'Ensembl/Havana gene',".
			"'name' => 'Merged Ensembl and Havana Genes',".
			"'label_key' => '[text_label] [display_label]',".
			"'default' => {'contigviewbottom' => 'transcript_label',".
			"'contigviewtop' => 'gene_label',".
			"'cytoview' => 'gene_label'},".
			"'key' => 'ensembl'}'";
		print "<ensembl> Updating web_data\n";
		update_analysis($caa, 'ensembl', undef, undef, $web_data);

	}

	if ($species =~ m/^(anopheles_gambiae)/) {

		print "<ensembl> Updating display_label to 'VectorBase gene'\n";
		update_analysis($caa, 'ensembl', undef, 'VectorBase gene');

		print "<anopheles_cdna_est> Switching displayable on\n";
		update_analysis($caa, 'anopheles_cdna_est', 1);

		print "<anopheles_cdna_est> Updating display_label to 'RNA (best)'\n";
		update_analysis($ofaa, 'anopheles_cdna_est', undef, 'RNA (best)');

	}
	
	if ($species =~ m/^(caenorhabditis_elegans)/) {

		# In C.elegans we have a couple of custom web_data columns that are 
		# perhaps most easily patched by reading from another logic_name 
		# - logic_name of 'ncRNA' has the same web_data as logic_name of 'tRNA'
		# - logic_name of 'Pseudogene' has the same web_data as logic_name of 'wormbase' 

		print "<$cdb> Overwriting web_data for logic_name 'ncRNA' ".
			"with web_data from 'tRNA'\n";
		my $tRNA = $caa->fetch_by_logic_name('tRNA');
		update_analysis($caa, 'ncRNA', undef, undef, $tRNA->web_data());

		print "<$cdb> Overwriting web_data for logic_name 'Pseudogene' ".
			"with web_data from 'wormbase'\n";
		my $wormbase = $caa->fetch_by_logic_name('wormbase');
		update_analysis($caa, 'Pseudogene', undef, undef, $wormbase->web_data());
		

	}
	

}

sub get_af_logic_names{

	my ($db, $molecule) = @_;

	my $sql = "select distinct logic_name from ".$molecule."_align_feature daf, analysis a ".
		"where daf.analysis_id=a.analysis_id;";
	return $db->dbc->db_handle->selectcol_arrayref($sql);

}

sub update_analysis {

	my ($aa, $logic_name, $displayable, $display_label, $web_data) = @_;

	my $analysis = $aa->fetch_by_logic_name($logic_name);
	throw("Analysis '$logic_name' is not defined") unless defined $analysis;

	if (defined $displayable) {
		print "\t[".$aa->db->dbc->dbname."] Updating '$logic_name' displayable from '".$analysis->displayable()."' to '".$displayable."'\n";
		$analysis->displayable($displayable)
	}
	if (defined $display_label) {
		print "\t[".$aa->db->dbc->dbname."] Updating '$logic_name' display_label from '".$analysis->display_label()."' to '".$display_label."'\n";
		$analysis->display_label($display_label);
	}
	if (defined $web_data) {
		print "\t[".$aa->db->dbc->dbname."] Updating '$logic_name' web_data from \"".$analysis->web_data()."\" to \"".$web_data."\"\n";
		$analysis->web_data($web_data);
	}

	$aa->update($analysis) if $update;

}

sub usage{
  exec('perldoc', $0);
  exit;
}

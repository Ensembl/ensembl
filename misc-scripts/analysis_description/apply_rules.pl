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
my $dbname = '';
my $help = 0;
my $version = 52;

&GetOptions (
	'host|dbhost=s'       => \$dbhost,
	'port|dbport=s'       => \$dbport,
	'user|dbuser=s'       => \$dbuser,
	'pass|dbpass=s'       => \$dbpass,
	'version=s'           => \$version,
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


$dsn = "DBI:mysql:database=" . $dbname . ";host=" . $dbhost . ";port=" . $dbport;
eval{ 
	$dbh = DBI->connect($dsn, $dbuser, $dbpass, 
						{'RaiseError' => 1,
						 'PrintError' => 0});
};


# get core databases;
my $sql = "show databases like '%core_$version%'";
my $cdbs  = $dbh->selectcol_arrayref($sql);
#print Dumper $cdbs;

foreach my $cdb (@$cdbs) {

	### implements rule 1. ###
	
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
				print("<$ln> exists only in otherfeatures, setting displayable 1 for otherfeatures\n");
				update_analysis($ofaa, $ln, 1);
	
			}
			delete $daf_logic_names{lc($ln)};

		}

		foreach my $ln (keys %daf_logic_names) {
				print("<$ln> exists only in core, setting displayable 1 for core\n");			
				update_analysis($caa, $ln, 1);

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
		$cdnaaa = $cdnadb->get_AnalysisAdaptor();

		my %alias = (
			'homo_sapiens' => 'Human',
			'mus_musculus' => 'Mouse'
		);

		my $ln = lc($alias{$species}).'_cdna';
		print "[$cdb] Switching off displayable for $ln\n";
		update_analysis($caa, $ln, 0);

		my $dl = $alias{$species}.' cDNA';
		print "[$cdnadb] Updating display_label for $ln to '$dl'\n";
		update_analysis($cdnaaa, $ln, 1, $dl);
			

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
		$vegaaa = $cdnadb->get_AnalysisAdaptor();
	
		my $vega_daf_logic_names = get_af_logic_names($vegadba, 'dna');
		my $vega_paf_logic_names = get_af_logic_names($vegadba, 'protein');
		
		foreach my $ln (@$vega_daf_logic_names, @$vega_paf_logic_names) {

			print "<$vegadb> Switching align_feature '$ln' displayable off\n";
			update_analysis($cdnaaa, $ln, 0);
			

		}

	}
	

}





sub get_af_logic_names{

	my ($db, $molecule) = @_;

	my $sql = "select distinct logic_name from ".$molecule."_align_feature daf, analysis a ".
		"where daf.analysis_id=a.analysis_id;";
	return $db->dbc->db_handle->selectcol_arrayref($sql);

}

sub update_analysis {

	my ($aa, $logic_name, $displayable, $display_label) = @_;

	my $analysis = $aa->fetch_by_logic_name($logic_name);

	if (defined $displayable) {
		print "\t[".$aa->db->dbc->dbname."] Updating '$logic_name' displayable from ".$analysis->displayable()." to ".$displayable."\n";
		$analysis->displayable($displayable)
	}
	if (defined $display_label) {
		print "\t[".$aa->db->dbc->dbname."] Updating '$logic_name' display_label from ".$analysis->display_label()." to ".$display_label."\n";
		$analysis->display_label($display_label);
	}

	$aa->update($analysis) if $update;

}

sub usage{
  exec('perldoc', $0);
  exit;
}

#!/opt/local/bin/perl
###!/usr/local/ensembl/bin/perl

# POD documentation - main docs before the code

=pod

=head1 NAME

  apply_rules.pl

=head1 SYNOPSIS

 script applies additonal rules to the analysis description table posterior to
 description loading. The script will create a file, debug_analysis_descriptions.txt
 with some verbose information about the different updates being done. This
 file should be kept for debugging purposes.

=head1 DESCRIPTION

 The rules (for details talk to Steve T.):
   1) For each species, dna_align_features and prediction_transcripts are displayed
      according to which database they can be retrieved from:
      - if they are present in both otherfeatures and core, then only the align_features
        from otherfeatures are to be displayed, ie we need to set the displayable 
        entry to 0 in core
      - prediction_transcripts are the opposite, ie off in other_features if they're on
        in core
      - if they are only present in one of these two databases then show them
        from that source

   2) Logic names of human_cdna and mouse_cdna are slightly different in that 
      they should be switched off (ie set the displayable entry to 0) in human 
      and mouse core databases respectively - they are superceded by the cDNA
      update features in the cDNA databases

   3) These cDNA_update features should have a display label of 'Mouse cDNA' in 
      the mouse_cdna database, and 'Human cDNA' in the human_cdna database. In addition
      they need another label (multi_caption) for multi location views

   4) All align_features from vega databases should be switched off

   5) Genes with a logic_name of 'otter' have a display label of 'Vega Havana gene' in
      Vega mouse and Vega human.

   6) Anopholes, human and mouse have different web_data columns from the definition file

   7) C.elegans has some unique web_data columns

   8) human and mouse have a different set of align_features switched on by default than
      is defined in the definition file (don't need the 'default'=>{'....'} entry

TO DO for e57:

- check that rat_cdna/rat_est are switched on in core and off in other_feature; if not add a rule
- add a rule to update web_data for GSTEN in other_features (see tetraodon_nigroviridis_otherfeatures_56_8b)
- add a rule / mdify analysis.descriptions for chimp otherfeatures to not display chimp_cdna, chimp_est 
  and human_cdna genes (see pan_troglodytes_otherfeatures_56_21l).

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

$| = 1;

my ($dsn, $dbh);
# Analysis adaptors
my ($caa, $ofaa, $cdnaaa, $vegaaa);

my $dbhost = '';
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;
my $help = 0;
my $version = 56;
my $update;
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

open (OUT,">debug_analysis_descriptions.txt") or throw("Failed to open file  debug_analysis_descriptions.txt");

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
    warn ("Displayable flag for analysis '$logic_name' has to be either 0 or 1, but not '$displayable'!")
      unless ($displayable =~ m/^[01]$/);
    throw("In the analysis_description file, logic name '$logic_name' should contain, at least, 5 columns: Number, logic_name, description, display_label and displayable. Fix it !!") unless (defined $displayable);

    #some entries in the file have no web_data defined, it is empty
    if (defined $web_data){
      #print OUT join("\t", $logic_name, $description, $display_label, $displayable, $web_data), "\n";
      $reference{lc($logic_name)} = {
	nr            => "$nr",
	description   => "$description",
	display_lable => "$display_label",
	displayable   => "$displayable",
	web_data      => "$web_data"
      }
    }
    else{
      #print OUT join("\t", $logic_name, $description, $display_label, $displayable), "\n";
      $reference{lc($logic_name)} = {
	nr            => "$nr",
	description   => "$description",
	display_lable => "$display_label",
	displayable   => "$displayable",
      }
    }
  }
  close FH;
	
} else {
  throw("Need to pass reference file with analysis descriptions!");
}


$dsn = "DBI:mysql:host=" . $dbhost . ";port=" . $dbport;
eval{ 
  $dbh = DBI->connect(
    $dsn, $dbuser, $dbpass, 
    {'RaiseError' => 1,
     'PrintError' => 0});
};


# get core database(s);
my $pat = defined $dbname ? $dbname : "%core_$version%";
my $sql = "show databases like '$pat'";
my $cdbs  = $dbh->selectcol_arrayref($sql);

foreach my $cdb (@$cdbs) {

  print OUT "\nStudying $cdb\n";
  (my $species = $cdb) =~ s/(.+)_core_${version}_\d+[a-z]$/$1/;

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
		
    print OUT "No otherfeatures db for " . $cdb . "! Setting all displayable entires to 1\n";
    my $daf_logic_names = get_af_logic_names($cdba, 'dna');
    map { update_analysis($caa, $_, 1) } @$daf_logic_names;
  }
  else {
    print OUT "Both core and otherfeatures dbs exist. Need to analyse dna_align_features ...\n";
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
	print OUT "<$ln> exists in both, setting displayable 0 for core and 1 for otherfeatures\n";
	update_analysis($caa, $ln, 0);
	update_analysis($ofaa, $ln, 1);	
      }
      else {	
	print OUT "<$ln> exists only in otherfeatures, setting displayable according to reference file\n";
	update_analysis($ofaa, $ln, $reference{lc($ln)}{displayable});	
      }
      delete $daf_logic_names{lc($ln)};
    }

    foreach my $ln (keys %daf_logic_names) {
      print OUT "<$ln> exists only in core, setting displayable according to reference file\n";
      update_analysis($caa, $ln, $reference{lc($ln)}{displayable});
    }

    # set prediction transcripts off in other_features if they're present in the core_db
    my %pt_logic_names;
    my $cpt_logic_names   = get_pt_logic_names($cdba);
    my $ofpt_logic_names  = get_pt_logic_names($ofdba);
    map {$pt_logic_names{lc($_)}++} (@$cpt_logic_names, @$ofpt_logic_names);
    foreach my $ln (keys %pt_logic_names) {
      if ($pt_logic_names{lc($ln)} == 2) {
	print OUT "<$ln> prediction transcript exists in both, setting displayable to 0 for otherfeatures\n";
	update_analysis($ofaa, $ln, 0);
      }
    }

    if ($species =~ m/^(homo_sapiens|mus_musculus)/) {

      ### new rule:in human and mouse, the only align_features on by default
      ### are cDNA update and CCDS
     
      my $core_daf_logic_names = get_af_logic_names($caa, 'dna');
      my $core_paf_logic_names = get_af_logic_names($caa, 'protein');
		
      foreach my $ln (@$core_daf_logic_names, @$core_paf_logic_names) {
	next if ($ln eq 'CCDS'); #this does not apply to CCDS
	my $ad = $caa->fetch_by_logic_name($ln);
	if ($ad->web_data() ne ''){
	  my $new_display_label = $ad->web_data();
	  delete $new_display_label->{'default'};
	  print OUT "<$cdb> Switching align_feature '$ln' from " . $caa->dump_data($ad->web_data()) . " to " . $caa->dump_data($new_display_label) ."\n";
	  update_analysis($caa, $ln, $ad->displayable(),$ad->display_label,$new_display_label);
	}
      }

      ### implements rule 2. and 3. ###

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
      print OUT "<$cdb> Switching off displayable for $ln\n";
      update_analysis($caa, $ln, 0);
      
      my $dl = $alias{$species}.' RefSeq/EMBL cDNA';
      my $ad = $caa->fetch_by_logic_name($ln);
      my $web_data = $ad->web_data();
      $web_data->{'multi_caption'} = 'Mouse/Human specific cDNA';
      $web_data->{'key'}  = 'species_specific_cdna';
      $web_data->{'name'} = $dl;
      print OUT "<$cdnadb> Updating display_label for cDNA_update to '$dl' and adding multi_caption and other entries to web_data\n";
      update_analysis($cdnaaa, 'cDNA_update', 1, $dl, $web_data);

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
	print OUT "<$vegadb> Switching align_feature '$ln' displayable off\n";
	update_analysis($vegaaa, $ln, 0);
      }
    
      ### implements rule 5. ###

      my $display_label = 'Vega Havana gene';
      print OUT "<$cdb> Updating display_label for logic_name otter to '$display_label'\n";
      update_analysis($vegaaa, 'otter', undef, $display_label);

      # 6) There are two more rules for the logic_name 'ensembl' in core databases 
      #    depending on the species:
      #	 - all but mouse, human and anopheles, use definition file.
      #	 - anopheles needs a display_label of 'VectorBase gene'
      #	 - mouse and human need a different web_data column:
     
      #web_data should be a hash reference
      $web_data = {
	'colour_key' => '[biotype]_[status]',
	'caption'    => 'Ensembl/Havana gene',
	'name'       => 'Merged Ensembl and Havana Genes',
	'label_key'  => '[text_label] [display_label]',
	'default'    => {
	  'contigviewbottom'     => 'transcript_label',
	  'contigviewtop'        => 'gene_label',
	  'cytoview'             => 'gene_label',
	  'MultiTop'             => 'gene_label',
	  'MultiBottom'          => 'collapsed_label',
	  'alignsliceviewbottom' => 'as_collapsed_label',
       },
	'multi_caption' => 'Ensembl genes, or Merged Ensembl and Havana genes',
	'key'           => 'ensembl',
      };
      print OUT "<ensembl> Updating web_data\n";
      update_analysis($caa, 'ensembl', undef, undef, $web_data);
    }

    if ($species =~ m/^(anopheles_gambiae)/) {
      
      print OUT "<ensembl> Updating display_label to 'VectorBase gene'\n";
      update_analysis($caa, 'ensembl', undef, 'VectorBase gene');

      print OUT "<anopheles_cdna_est> Switching displayable on\n";
      update_analysis($caa, 'anopheles_cdna_est', 1);

      print OUT "<anopheles_cdna_est> Updating display_label to 'RNA (best)'\n";
      update_analysis($ofaa, 'anopheles_cdna_est', undef, 'RNA (best)');
    }
	
    #rule 7
    if ($species =~ m/^(caenorhabditis_elegans)/) {
      
      # In C.elegans we have a couple of custom web_data columns that are 
      # perhaps most easily patched by reading from another logic_name 
      # - logic_name of 'ncRNA' has the same web_data as logic_name of 'tRNA'
      # - logic_name of 'Pseudogene' has the same web_data as logic_name of 'wormbase' 
      
      print OUT "<$cdb> Overwriting web_data for logic_name 'ncRNA' ".
	"with web_data from 'tRNA'\n";
      my $tRNA = $caa->fetch_by_logic_name('tRNA');
      update_analysis($caa, 'ncRNA', undef, undef, $tRNA->web_data());
      
      print OUT "<$cdb> Overwriting web_data for logic_name 'Pseudogene' ".
	"with web_data from 'wormbase'\n";
      my $wormbase = $caa->fetch_by_logic_name('wormbase');
      update_analysis($caa, 'Pseudogene', undef, undef, $wormbase->web_data());
    }
  }
}
close OUT;

sub get_af_logic_names{
  my ($db, $molecule) = @_;
  my $sql = "select distinct logic_name from ".$molecule."_align_feature daf, analysis a ".
    "where daf.analysis_id=a.analysis_id;";
  return $db->dbc->db_handle->selectcol_arrayref($sql);
}

sub get_pt_logic_names{
  my ($db) = @_;
  my $sql = "select distinct logic_name from prediction_transcript pt, analysis a where pt.analysis_id=a.analysis_id;";
  return $db->dbc->db_handle->selectcol_arrayref($sql);
}

sub update_analysis {
  my ($aa, $logic_name, $displayable, $display_label, $web_data) = @_;
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  throw("Analysis '$logic_name' is not defined") unless defined $analysis;
  if (defined $displayable) {
    print OUT "\t[".$aa->db->dbc->dbname."] Updating '$logic_name' displayable from '".$analysis->displayable()."' to '".$displayable."'\n";
    $analysis->displayable($displayable)
  }
  if (defined $display_label) {
    print OUT "\t[".$aa->db->dbc->dbname."] Updating '$logic_name' display_label from '".$analysis->display_label()."' to '".$display_label."'\n";
    $analysis->display_label($display_label);
  }
  if (defined $web_data) {
    print OUT "\t[".$aa->db->dbc->dbname."] Updating '$logic_name' web_data from \"".$aa->dump_data($analysis->web_data())."\" to \"".$aa->dump_data($web_data)."\"\n";
    $analysis->web_data($web_data);
  }
  $aa->update($analysis) if $update;
}

sub usage{
  exec('perldoc', $0);
  exit;
}

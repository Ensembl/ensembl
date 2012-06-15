use strict;

# Sets display_xref_ids for novel genes in the "to" database based
# on their orthologs in the "from" database. Can also project GO xrefs.
# Orthology relationships are read from a Compara database.

use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::OntologyTermAdaptor;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end eprof_dump);

my $method_link_type = "ENSEMBL_ORTHOLOGUES";

my %seen;

my ($conf, $registryconf, $version, $compara, $from_species, @to_multi, $print, $names, $go_terms, $delete_names, $delete_go_terms, $no_backup, $full_stats, $descriptions, $release, $no_database, $quiet, $max_genes, $one_to_many, $go_check, $all_sources, $delete_only,  $to_species, $from_gene);

GetOptions('conf=s'          => \$conf,
	   'registryconf=s'  => \$registryconf,
	   'version=i'       => \$version,
	   'compara=s'       => \$compara,
	   'from=s'          => \$from_species,
	   'to=s'            => \@to_multi,
	   'method=s'        => \$method_link_type,
	   'names'           => \$names,
	   'go_terms'        => \$go_terms,
	   'print'           => \$print,
	   'delete_names'    => \$delete_names,
	   'delete_go_terms' => \$delete_go_terms,
	   'nobackup'        => \$no_backup,
	   'full_stats'      => \$full_stats,
       'descriptions'    => \$descriptions,
	   'release=i'       => \$release,
	   'no_database'     => \$no_database,
	   'quiet'           => \$quiet,
	   'max_genes=i'     => \$max_genes,
	   'one_to_many'     => \$one_to_many,
	   'go_check'        => \$go_check,
	   'all_sources'     => \$all_sources,
	   'delete_only'     => \$delete_only,
	   'help'            => sub { usage(); exit(0); });

$| = 1; # auto flush stdout

$descriptions = 1;

if (!$conf && !$registryconf && !$delete_only) {

  print STDERR "Configuration file must be supplied via -conf or -registryconf argument\n";
  usage();
  exit(1);

} elsif (!$from_species && !$delete_only)  {

  print STDERR "From species must be supplied via -from argument\n";
 usage();
  exit(1);

} elsif (!@to_multi) {

  print STDERR "At least one target species must be supplied via the -to argument\n";
  usage();

  exit(1);

} elsif (!$release && !$no_database) {

  print STDERR "Release must be specified via -release argument unless -no_database is used\n";
  usage();
  exit(1);

}

if (!$go_terms && !$names && !$delete_only) {

  print STDERR "One or both of --names or --go_terms must be specified unless only -delete_only is being used\n";
  print STDERR "Use --help for more detailed usage informaion\n";
  exit(1);

}

# only certain types of homology are considered
my @homology_types_allowed = ("ortholog_one2one","apparent_ortholog_one2one");
if ($one_to_many) {
  push @homology_types_allowed, "ortholog_one2many","apparent_ortholog_one2many";
}

# only these evidence codes will be considered for GO term projection
my @evidence_codes = ( "IDA", "IEP", "IGI", "IMP", "IPI", "EXP" );

#  IC Inferred by curator
#  IDA Inferred from direct assay
#  IEA Inferred from electronic annotation
#  IGI Inferred from genetic interaction
#  IMP Inferred from mutant phenotype
#  IPI Inferred from physical interaction
#  ISS Inferred from sequence or structural similarity
#  NAS Non-traceable author statement
#  ND No biological data available
#  RCA Reviewed computational analysis
#  TAS Traceable author statement

@to_multi = split(/,/,join(',',@to_multi));

# load from database and conf file
Bio::EnsEMBL::Registry->no_version_check(1);


# Registryconf is either the registry configuration passed from the submit_projections.pl 
# script or a file name containing the same information that is passed on the command line.

my $args;

if (defined($registryconf)) {
	if (-f $registryconf) {
		open(CONF, $registryconf);
		my @contents = <CONF>;
		$args = eval(join("\n", @contents));
		close(CONF);
	} else {
		$args = eval($registryconf);
	}
}

Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs(@{$args});

Bio::EnsEMBL::Registry->load_all($conf, 0, 1); # options mean "not verbose" and "don't clear registry"

# only delete names/GO terms if -delete_only has been specified
if ($delete_only) {

    print "Just deleting, no projection\n";
    foreach my $to_species (@to_multi) {

    	my $to_ga  = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'Gene');
    	die("Can't get gene adaptor for $to_species - check database connection details; make sure meta table contains the correct species alias\n") if (!$to_ga);
    	delete_names($to_ga) if ($delete_names);
    	delete_go_terms($to_ga) if ($delete_go_terms);
    }

    exit(0);
}

# Get Compara adaptors - use the one specified on the command line, or the first one
# defined in the registry file if not specified

my $mlssa;
my $ha;
my $ma;
my $gdba;

if ($compara) {

   $mlssa = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'MethodLinkSpeciesSet');
   $ha    = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'Homology');
   $ma    = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'Member');
   $gdba  = Bio::EnsEMBL::Registry->get_adaptor($compara, "compara", "GenomeDB");

   die "Can't connect to Compara database specified by $compara - check command-line and registry file settings" if (!$mlssa || !$ha || !$ma ||!$gdba);

} else {

   $mlssa = @{Bio::EnsEMBL::Registry->get_all_adaptors(-group => "compara", -type => "MethodLinkSpeciesSet")}[0];
   $ha    = @{Bio::EnsEMBL::Registry->get_all_adaptors(-group => "compara", -type => "Homology")}[0];
   $ma    = @{Bio::EnsEMBL::Registry->get_all_adaptors(-group => "compara", -type => "Member")}[0];
   $gdba  = @{Bio::EnsEMBL::Registry->get_all_adaptors(-group => "compara", -type => "GenomeDB")}[0];

   die "Can't connect to Compara database from registry - check registry file settings" if (!$mlssa || !$ha || !$ma ||!$gdba);

}



my $from_ga = Bio::EnsEMBL::Registry->get_adaptor($from_species, 'core', 'Gene');

my %projections_by_evidence_type;
my %projections_by_source;
my %forbidden_terms;

my $from_mammal;
my $to_mammal;

foreach my $local_to_species (@to_multi) {

  $to_species = $local_to_species;
  my $to_ga   = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'Gene');
  die("Can't get gene adaptor for $to_species - check database connection details; make sure meta table contains the correct species alias\n") if (!$to_ga);
  my $to_dbea = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'DBEntry');


    # Determine if source and target species are both mammals or not: class Mammalia in meta
    # table, meaning full GO projection is allowable
    # Platypus is not treated as a mammal here, but "is_therian" doesn't mean anything.
    
    $from_mammal = is_mammal($from_species);
    $to_mammal = is_mammal($to_species);
    
    if ($from_mammal && !$to_mammal) {
        # Determine GO terms that are not applicable for the present non-mammalian species.
        # Consult ontology for list of descendent terms.
        # Using a temporary list of general forbidden terms. Presently consists of four GO terms
        # recommended by GOA. e.g. "walks on legs" should not be projected to fish
        %forbidden_terms = get_ontology_terms('GO:0032501','GO:0007610','GO:0048856','GO:0051704');
        
        # The forbidden_terms list is used in unwanted_go_term();
    }
    
    #######

  write_to_projection_db($to_ga->dbc(), $release, $from_species, $from_ga->dbc(), $to_species) unless ($no_database);

  backup($to_ga) if (!$no_backup);

  delete_names($to_ga) if ($delete_names);
  delete_go_terms($to_ga) if ($delete_go_terms);

  # get taxonomy ids from core databases
  my $meta_container =  Bio::EnsEMBL::Registry->get_adaptor($from_species, 'Core', 'MetaContainer');
  my $from_taxon_id = $meta_container->get_taxonomy_id();
  $meta_container =  Bio::EnsEMBL::Registry->get_adaptor($to_species, 'Core', 'MetaContainer');
  my $to_taxon_id = $meta_container->get_taxonomy_id();

  # build Compara GenomeDB objects
  my $from_GenomeDB = $gdba->fetch_by_taxon_id($from_taxon_id);
  my $to_GenomeDB = $gdba->fetch_by_taxon_id($to_taxon_id);

  my $mlss = $mlssa->fetch_by_method_link_type_GenomeDBs($method_link_type, [$from_GenomeDB, $to_GenomeDB]);

  # get homologies from compara - comes back as a hash of arrays
  my $homologies = fetch_homologies($ha, $mlss, $from_species);

  my $str = "gene display_xrefs and descriptions" if ($names);
  $str .= " and " if ($names && $go_terms);
  $str .= "GO terms" if ($go_terms);

  print "Projecting $str from $from_species to $to_species\n";

  # build hash of external db name -> ensembl object type mappings
  my %db_to_type = build_db_to_type($to_ga);

  print "$to_species, before projection: \n";
  print_stats($to_ga);

  my $i = 0;
  my $total_genes = scalar(keys %$homologies);
  my $last_pc = -1;

  print "Percentage complete: " if (!$quiet);

  foreach my $from_stable_id (keys %$homologies) {

    $last_pc = print_progress($i, $total_genes, $last_pc) if (!$quiet);

    last if ($max_genes && $i > $max_genes);

    $i++;

    $from_gene = $from_ga->fetch_by_stable_id($from_stable_id);


    next if (!$from_gene);
    # check for uniprot IDs in canonical transcript
    #TODO: This feature is waiting to be implemented!

    my @to_genes = @{$homologies->{$from_stable_id}};
    my $i = 1;
    foreach my $to_stable_id (@to_genes) {

      my $to_gene = $to_ga->fetch_by_stable_id($to_stable_id);

      next if (!$to_gene);

      project_display_names($to_ga, $to_dbea, $from_gene, $to_gene, $i, scalar(@to_genes), %db_to_type) if ($names);

      project_go_terms($to_ga, $to_dbea, $ma, $from_gene, $to_gene) if ($go_terms);

      $i++;

    }

  }

  if ($go_terms) {
    print "\nCleaning up GO projections ...";
    clean_up($to_ga->dbc());
  }

  print "\n$to_species, after projection: \n";
  print_stats($to_ga);

  print_full_stats() if ($full_stats);

}

# mammal here is a pseudonym for normal mammals and marsupials but not platypus.
sub is_mammal {
    my $species = shift;
    my $meta_container = Bio::EnsEMBL::Registry->get_adaptor($species,'core','MetaContainer');
    my @classifications = @{ $meta_container->list_value_by_key('species.classification') };

    foreach my $rank (@classifications) {
        if ($rank eq "Theria") {return 1;} 
    }
    return;
}

sub get_ontology_terms {
    my @starter_terms = @_;
    my $ontology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','Ontology','OntologyTerm');
    my %terms;
    foreach my $text_term (@starter_terms) {
        my $ont_term = $ontology_adaptor->fetch_by_accession($text_term);
        my $term_list = $ontology_adaptor->fetch_all_by_ancestor_term($ont_term);
        foreach my $term (@{$term_list}) {
            $terms{$term->accession} = 1;
        }
    }
    return %terms;
}

# --------------------------------------------------------------------------------

sub project_display_names {

  my ($to_ga, $to_dbea, $from_gene, $to_gene,  $gene_number, $total_gene_number, %db_to_type) = @_;

  my $dbEntry = $from_gene->display_xref();
  my $to_source = $to_gene->display_xref()->dbname() if ($to_gene->display_xref());
  my $from_source = $from_gene->display_xref()->dbname() if ($from_gene->display_xref());

  my $from_latin_species = ucfirst(Bio::EnsEMBL::Registry->get_alias($from_species));

  # if no display name set, do the projection
  if (check_overwrite_display_xref($to_gene, $from_source, $to_source, $dbEntry, $gene_number)) {

    if ($dbEntry) {

      my $dbname = $dbEntry->dbname();

      return if (!$all_sources && $dbname !~ /HGNC/);

      # Skip clone names if projecting all sources
      return if (lc($dbname) =~ /clone/);

      # Modify the dbEntry to indicate it's not from this species - set info_type & info_text
      my $info_txt = "from $from_latin_species gene " . $from_gene->stable_id();

      # modify the display_id to have " (1 of 3)" etc if this is a one-to-many ortholog
      my $tuple_txt = "";
      if ($total_gene_number > 1) {
          $tuple_txt = " ($gene_number of $total_gene_number)";
          my $existing = $dbEntry->display_id();
          $existing =~ s/ \(\d+ of \d+\)//;
          $dbEntry->display_id($existing . $tuple_txt);
          $info_txt .= $tuple_txt;
      }

      # Add description to the gene if required
      # This should probably be in another column in the xref table but we just use the gene
      # description column for now.
      $to_gene->description($from_gene->description()) if ($descriptions && $from_gene->description());

      $dbEntry->info_type("PROJECTION");
      $dbEntry->info_text($info_txt);

      # Add the xref to the "to" gene, or transcript or translation depending on what the
      # other xrefs from this dbname as assigned to (see build_db_to_type)
      # Note that if type is not found, it means that we're dealing with a db that has no
      # xrefs in the target database, e.g. MarkerSymbol in mouse -> rat
      # In this case just assign to genes

      my @to_transcripts = @{$to_gene->get_all_Transcripts};
      my $to_transcript = $to_transcripts[0];
      
      # Force loading of external synonyms for the xref
      $dbEntry->get_all_synonyms();

      $dbname = $dbEntry->dbname();

      my $type = $db_to_type{$dbname};

      if ($type eq "Gene" || $dbname eq 'HGNC' || !$type) {

	$to_gene->add_DBEntry($dbEntry);
	$to_dbea->store($dbEntry, $to_gene->dbID(), 'Gene', 1) if (!$print);

      } elsif ($type eq "Transcript" || $dbname eq 'HGNC_transcript_name') {
	
	$to_transcript->add_DBEntry($dbEntry);
	$to_dbea->store($dbEntry, $to_transcript->dbID(), 'Transcript', 1) if (!$print);

      } elsif ($type eq "Translation") {

	my $to_translation = $to_transcript->translation();
	return if (!$to_translation);
	$to_translation->add_DBEntry($dbEntry);
	$to_dbea->store($dbEntry, $to_translation->dbID(), 'Translation',1) if (!$print);

      } else {

	warn("Can't deal with xrefs assigned to $type (dbname=" . $dbEntry->dbname . ")\n");
	return;

      }

      # Set gene status to "KNOWN_BY_PROJECTION" and update display_xref
      # also set the status of the gene's transcripts
      $to_gene->status("KNOWN_BY_PROJECTION");
      $to_gene->display_xref($dbEntry);
      foreach my $transcript (@{$to_gene->get_all_Transcripts()}) {
	$transcript->status("KNOWN_BY_PROJECTION");
      }
	
      print $to_gene->stable_id() . " --> " . $dbEntry->display_id() . "\n" if ($print);

      # update the gene so that the display_xref_id is set
      $to_ga->update($to_gene) if (!$print);

      # keep track of where each projection came from
      $projections_by_source{$dbEntry->dbname()}++;

    }

  }

}

# --------------------------------------------------------------------------------

sub project_go_terms {

  my ($to_ga, $to_dbea, $ma, $from_gene, $to_gene) = @_;
# TODO: compara member adaptor $ma doesn't appear to be used in this method

  # GO xrefs are linked to translations, not genes
  # Project GO terms between the translations of the canonical transcripts of each gene
  my $from_translation = get_canonical_translation($from_gene);
  my $to_translation   = get_canonical_translation($to_gene);

  return if (!$from_translation || !$to_translation);

  my $from_latin_species = ucfirst(Bio::EnsEMBL::Registry->get_alias($from_species));

  my $to_go_xrefs = $to_translation->get_all_DBEntries("GO") if ($go_check);

 DBENTRY: foreach my $dbEntry (@{$from_translation->get_all_DBEntries("GO")}) { 

    next if (!$dbEntry || $dbEntry->dbname() ne "GO" || ref($dbEntry) ne "Bio::EnsEMBL::OntologyXref");

    # Skip the whole dbEntry if one or more if its evidence codes isn't in the whitelist
    foreach my $et (@{$dbEntry->get_all_linkage_types}){
      next DBENTRY if (!grep(/$et/, @evidence_codes));
    }
    next if (unwanted_go_term($dbEntry));

    # check that each from GO term isn't already projected
    next if ($go_check && go_xref_exists($dbEntry, $to_go_xrefs));
    
    # Force loading of external synonyms for the xref
    $dbEntry->get_all_synonyms();

    # record statistics by evidence type
    foreach my $et (@{$dbEntry->get_all_linkage_types}){
      $projections_by_evidence_type{$et}++;
    }

    # Change linkage_type for projection to IEA (in the absence of a specific one for projections)
    $dbEntry->flush_linkage_types();
    $dbEntry->add_linkage_type("IEA");

    my $txt = "from $from_latin_species translation " . $from_translation->stable_id();
    $dbEntry->info_type("PROJECTION");
    $dbEntry->info_text($txt);

    $to_translation->add_DBEntry($dbEntry);

    print $from_gene->stable_id() . " " . $from_translation->stable_id() . " " .  $dbEntry->display_id() . " --> " . $to_gene->stable_id() . " " . $to_translation->stable_id() . "\n" if ($print);

    $to_dbea->store($dbEntry, $to_translation->dbID(), 'Translation', 1) if (!$print);

  }

}

sub unwanted_go_term {
    my $dbEntry = shift;
    
    if ($from_mammal && !$to_mammal) {
        if (exists $forbidden_terms{$dbEntry->primary_id}) {
            return 1;
        }
    }
    return 0;
}

# --------------------------------------------------------------------------------

sub go_xref_exists {

  my ($dbEntry, $to_go_xrefs) = @_;

  foreach my $xref (@{$to_go_xrefs}) {

    next if (ref($dbEntry) ne "Bio::EnsEMBL::OntologyXref" || ref($xref) ne "Bio::EnsEMBL::OntologyXref");

    if ($xref->dbname() eq $dbEntry->dbname() &&
	$xref->primary_id() eq $dbEntry->primary_id() &&
	join("", @{$xref->get_all_linkage_types()}) eq join("", @{$dbEntry->get_all_linkage_types()})) {
      return 1;
    }

    # if a GO term with the same accession, but IEA evidence code, exists, also don't project, as this
    # will lead to duplicates when the projected term has its evidence code changed to IEA after projection
    if ($xref->primary_id() eq $dbEntry->primary_id()) {
      foreach my $evidence_code (@{$xref->get_all_linkage_types()}) {
	return 1 if ($evidence_code eq "IEA");
      }
    }

  }

  return 0;

}

# ----------------------------------------------------------------------

sub print_stats {

  my ($to_ga) = @_;

  my $total_genes = count_rows($to_ga, "SELECT COUNT(*) FROM gene g");

  my $count;

  if ($names) {

    $count = count_rows($to_ga, "SELECT COUNT(*) FROM gene g, xref x WHERE g.display_xref_id=x.xref_id AND g.display_xref_id IS NOT NULL AND  (x.info_type != 'PROJECTION' || x.info_type IS NULL)");
    printf("Gene names: unprojected %d (%3.1f\%)" , $count, (100 * $count / $total_genes));

    my $projected = count_rows($to_ga, "SELECT COUNT(*) FROM gene g, xref x WHERE g.display_xref_id=x.xref_id AND x.info_type='PROJECTION'");
    printf(" projected %d (%3.1f\%)" , $projected, (100 * $projected / $total_genes));

    $count = count_rows($to_ga, "SELECT COUNT(*) FROM gene g, xref x, external_db e WHERE g.display_xref_id=x.xref_id AND x.external_db_id=e.external_db_id AND e.db_name IN ('RefSeq_mRNA_predicted', 'RefSeq_ncRNA_predicted', 'RefSeq_peptide_predicted')");
    printf(" predicted %d (%3.1f\%)" , $count, (100 * $count / $total_genes));

    $count = count_rows($to_ga, "SELECT COUNT(*) FROM gene g WHERE display_xref_id IS NOT NULL");
    printf(" total genes with names %d (%3.1f\%)\n" , $count, (100 * $count / $total_genes));

    if ($projected > 0) {
      my $one2many = count_rows($to_ga, "SELECT COUNT(*) FROM gene g, xref x WHERE g.display_xref_id=x.xref_id AND x.info_type='PROJECTION' AND x.display_label LIKE '%(% of %)%'");
      my $one2one = $projected - $one2many;
      printf("Of the %d projected genes, %d (%3.1f\%) are from one-one mappings, %d (%3.1f\%) from one-many mappings\n", $projected, $one2one, (100 * $one2one/$projected), $one2many, (100 * $one2many / $projected));
    }
  }

  if ($go_terms) {

    print "GO xrefs: total ";
    print &count_rows($to_ga, "SELECT COUNT(DISTINCT(x.dbprimary_acc)) FROM xref x, external_db e WHERE e.external_db_id=x.external_db_id AND e.db_name='GO'");

    print " projected ";
    print &count_rows($to_ga, "SELECT COUNT(DISTINCT(x.dbprimary_acc)) FROM xref x, external_db e WHERE e.external_db_id=x.external_db_id AND e.db_name='GO' AND x.info_type='PROJECTION'");

    print "\n";

  }

}

# ----------------------------------------------------------------------

sub print_full_stats {

# GO terms
  if ($go_terms) {

    print "\nProjected terms by evidence code:\n";
    my $total;
    foreach my $et (sort keys %projections_by_evidence_type) {
	
      next if (!grep(/$et/, @evidence_codes));
	
      if ($et) {
	print $et . "\t" . $projections_by_evidence_type{$et} . "\n";
	$total += $projections_by_evidence_type{$et};
      }
    }
    print "Total:\t$total\n";
  }

  # display names
  if ($names) {

    print "\nProjected display names by source:\n";

    my $total;
    foreach my $source (sort keys %projections_by_source) {

      $total += $projections_by_source{$source};

    }

    foreach my $source (sort keys %projections_by_source) {

      my $n = $projections_by_source{$source};
      printf ("%s\t%d (%3.1f\%)\n",  $source, $n, (100 * $n / $total));

  }

  print "Total:\t$total\n";

  }


}

# ----------------------------------------------------------------------

sub count_rows {

  my ($adaptor, $sql) = @_;

  my $sth = $adaptor->dbc->prepare($sql);
  $sth->execute();

  return ($sth->fetchrow_array())[0];

}

# ----------------------------------------------------------------------

# create a hash of external_db_name -> ensembl_object_type
# used to assign projected xrefs to the "correct" type

sub build_db_to_type {

  my ($to_ga) = @_;

  my %db_to_type = ();

  my $sth = $to_ga->dbc()->prepare("SELECT DISTINCT e.db_name, ox.ensembl_object_type FROM external_db e, xref x, object_xref ox WHERE x.xref_id=ox.xref_id AND e.external_db_id=x.external_db_id");
  $sth->execute();
  my ($db_name, $type);
  $sth->bind_columns(\$db_name, \$type);
  while($sth->fetch()){
    $db_to_type{$db_name} = $type;
  }
  $sth->finish;

  return %db_to_type;

}

# ----------------------------------------------------------------------

sub delete_names {

  my ($to_ga) = @_;

  print "Setting gene display_xrefs and descriptions that were projected to NULL, and status to NOVEL\n";
  my $sth = $to_ga->dbc()->prepare("UPDATE gene g, xref x SET g.display_xref_id = NULL, g.description=NULL, g.STATUS='NOVEL' WHERE g.display_xref_id=x.xref_id AND x.info_type='PROJECTION'");
  $sth->execute();

  print "Deleting projected xrefs, object_xrefs and synonyms\n";
  $sth = $to_ga->dbc()->prepare("DELETE es FROM xref x, external_synonym es WHERE x.xref_id=es.xref_id AND x.info_type='PROJECTION'");
  $sth->execute();
  # avoid deleting projected GO terms - only want to delete the names here
  $sth = $to_ga->dbc()->prepare("DELETE x, ox FROM xref x, object_xref ox, external_db e WHERE x.xref_id=ox.xref_id AND x.external_db_id=e.external_db_id AND x.info_type='PROJECTION' AND e.db_name!='GO'");
  $sth->execute();

}

# ----------------------------------------------------------------------

sub delete_go_terms {

  my ($to_ga) = @_;

  print "Deleting projected GO terms\n";

  my $sth = $to_ga->dbc()->prepare("DELETE x, ox, gx FROM xref x, external_db e, object_xref ox, ontology_xref gx WHERE x.xref_id=ox.xref_id AND x.external_db_id=e.external_db_id AND ox.object_xref_id=gx.object_xref_id AND e.db_name='GO' AND x.info_type='PROJECTION'");
  $sth->execute();

  # note don't need to delete synonyms as GO terms don't have any
  # Also no effect on descriptions or status
}

# ----------------------------------------------------------------------

# Decide if a gene name should be overwritten
# Criteria: overwrite if:
#    - no existing display_xref
# or
#    - existing display_xref is RefSeq_*_predicted
#      AND from_gene is from "best" source external db,
#      e.g. HGNC in human, MGI in mouse

sub check_overwrite_display_xref {

  my ($to_gene, $from_dbname, $to_dbname, $ref_dbEntry, $gene_number) = @_;

  return 1 if (!$to_gene->external_name() && $to_species ne "zebrafish");

  if ($to_dbname eq "RefSeq_mRNA_predicted" || $to_dbname eq "RefSeq_ncRNA_predicted" || $to_dbname eq "RefSeq_peptide_predicted") {

    if (($from_species eq "human" && $from_dbname =~ /HGNC/) ||
    ($from_species eq "mouse" && $from_dbname =~ /MarkerSymbol/)) {

    if ($to_species eq "zebrafish" and is_in_blacklist($from_gene->display_xref)){
        return 0;
    }
      return 1;

    }

  }
  elsif ($to_species eq "zebrafish"){

    my $to_dbEntry = $to_gene->display_xref();
    my $from_dbEntry = $from_gene->display_xref();
    my $to_seq_region_name = $to_gene->seq_region_name();

    return 1 if ($to_dbname eq "Clone_based_ensembl_gene" or $to_dbname eq "Clone_based_vega_gene");

    my $name = $ref_dbEntry->display_id;
    $name =~ /(\w+)/; # remove (x of y) in name.
    $name = $1;

    if(defined($seen{$name})){
        $name = $seen{$name};
    }
    if ( $name =~ /C(\w+)orf(\w+)/){
        my $new_name = "C".$to_seq_region_name."H".$1."orf".$2;
        $seen{$new_name} = $name;
        $ref_dbEntry->display_id($new_name);
        return 1;
    }

    if (!defined ($to_dbEntry) || (($to_dbEntry->display_id =~ /:/) and $to_dbname eq "ZFIN_ID") ){
      if (is_in_blacklist($from_dbEntry)){
	      return 0;
      }
      else{
	      return 1;
      }
    }
  }
  return 0;

}


sub is_in_blacklist{
    #catches clones and analyses when projecting display xrefs.
    my ($dbentry) = shift;

    if (($dbentry->display_id =~ /KIAA/) || ( $dbentry->display_id =~ /LOC/)){
       return 1; # return yes that have found gene names that match the regular expression
    }
    elsif ($dbentry->display_id =~ /\-/){
       return 1;
    }
    elsif ($dbentry->display_id =~ /\D{2}\d{6}\.\d+/){
       #print "black listed item found ".$dbentry->display_id."\n";
        return 1;
    }
    else{
        return 0;
    }
    
}
# ----------------------------------------------------------------------

sub backup {

  my ($to_ga) = @_;

  my $dbc = $to_ga->dbc();
  my $host = $dbc->host();
  my $port = $dbc->port();
  my $user = $dbc->username();
  my $pass = $dbc->password();
  my $dbname = $dbc->dbname();

  foreach my $table ("gene", "xref", "object_xref") {
    unless (system("mysql -h$host -P$port -u$user -p$pass -N -e 'select * from $table' $dbname > $dbname.$table.backup; gzip -6 -f $dbname.$table.backup") == 0) {
      print STDERR "Can't dump the original $table table from $dbname for backup\n";
      exit 1;
    } else {
      print "Original $table table backed up in $dbname.$table.backup\n";
    }
  }

}

# ----------------------------------------------------------------------

sub write_to_projection_db {

  my ($to_dbc, $release, $from_species, $from_dbc, $to_species) = @_;

  my $host = $to_dbc->host();
  my $port = $to_dbc->port();
  my $user = $to_dbc->username();
  my $pass = $to_dbc->password();
  my $dbname = "projection_info";

  my $from_dbname = $from_dbc->dbname();
  my $to_dbname = $to_dbc->dbname();

  my $db = DBI->connect("dbi:mysql:host=$host;port=$port;database=$dbname", "$user", "$pass", {'RaiseError' => 1}) || die "Can't connect to " . $dbname;

  my $sth = $db->prepare("INSERT INTO projections (db_release, timestamp, from_db, from_species_latin, from_species_common, to_db, to_species_latin, to_species_common) VALUES (?, NOW(), ?, ?, ?, ?, ?, ?)");

  $sth->execute($release,
		$from_dbname,
		Bio::EnsEMBL::Registry->get_alias($from_species),
		$from_species,
		$to_dbname,
		Bio::EnsEMBL::Registry->get_alias($to_species),
		$to_species) || die "Can't write to projection info database\n";

  $sth->finish();

}

# ----------------------------------------------------------------------
# Fetch the homologies from the Compara database.
# Returns a hash of arrays:
# Key = "from" stable ID, value = array of "to" stable IDs

sub fetch_homologies {

  my ($ha, $mlss, $from_species) = @_;

  print "Fetching Compara homologies\n";

  my $from_species_alias = lc(Bio::EnsEMBL::Registry->get_alias($from_species));

  my %homology_cache;

  my $count = 0;

  my $homologies = $ha->fetch_all_by_MethodLinkSpeciesSet($mlss);

  foreach my $homology (@{$homologies}) {

    next if (!homology_type_allowed($homology->description));

    my @mas = @{$homology->get_all_Member_Attribute};

    # order of member-attributes is arbitrary, so need to find which one corresponds to the "from" species
    my @to_stable_ids;
    my $from_stable_id;

    foreach my $ma (@mas) {

      my ($member, $attribute) = @{$ma};

      if (lc($member->genome_db()->name()) eq $from_species_alias) {
          $from_stable_id = $member->stable_id();
      } else {
          push @to_stable_ids, $member->stable_id();
      }
    }

    print "Warning: can't find stable ID corresponding to 'from' species ($from_species_alias)\n" if (!$from_stable_id);

    push @{$homology_cache{$from_stable_id}}, @to_stable_ids;

    $count++;

  }

  print "Fetched " . $count . " homologies\n";

  return \%homology_cache;
}

# ----------------------------------------------------------------------

sub homology_type_allowed {

  my $h = shift;

  foreach my $allowed (@homology_types_allowed) {
    return 1 if ($h eq $allowed);
  }

  return undef;

}

# ----------------------------------------------------------------------
# Get the translation associated with the gene's canonical transcript

sub get_canonical_translation {

  my $gene = shift;

  my $canonical_transcript = $gene->canonical_transcript();

  if (!$canonical_transcript) {
    warn("Can't get canonical transcript for " . $gene->stable_id() . ", skipping this homology");
    return undef;
  }

  return $canonical_transcript->translation();;

}

# ----------------------------------------------------------------------

sub print_progress {

  my ($i, $total, $last_pc) = @_;

  my $pc = int ((100 * $i) / $total);

  if ($pc > $last_pc) {
    print "$pc ";
    $last_pc = $pc;
  }

  return $last_pc;

}

# ----------------------------------------------------------------------

sub clean_up {

  my ($to_dbc) = @_;

  $to_dbc->do("DROP TABLE IF EXISTS tmp_gx");

  $to_dbc->do("CREATE TEMPORARY TABLE tmp_gx SELECT gx.object_xref_id FROM ontology_xref gx LEFT JOIN object_xref ox ON ox.object_xref_id=gx.object_xref_id WHERE ox.object_xref_id IS NULL");

  $to_dbc->do("DELETE FROM ontology_xref WHERE object_xref_id IN (SELECT object_xref_id FROM tmp_gx)");

  $to_dbc->do("DROP TABLE tmp_gx");

}

# ----------------------------------------------------------------------



sub usage {

  print << "EOF";

  Sets display_xref_ids and/or GO terms for novel genes in the "to" database
  based on their orthologs in the "from" database. Orthology relationships
  are read from a Compara database.

 perl project_display_xrefs.pl {options}

 Options ([..] indicates optional):

  [--conf filepath]     the Bio::EnsEMBL::Registry configuration file. If none
                        is given, the one set in ENSEMBL_REGISTRY will be used
                        if defined, if not ~/.ensembl_init will be used.
                        Note only the Compara database needs to be defined here,
                        assuming the rest of the databases are on the server
                        defined by --registryconf
   

   --registryconf       There are two ways in which the registry configuration
                        information can be passed to the script. This information
		        is a hash that encodes the registry configuration parameters
			and can be passed as a string in a file or as a string on the 
			commandline. 


   --version

                        Note that a combination of the host/user and conf files
                        can be used. Databases specified in both will use the
                        conf file setting preferentially.

   --from string        The species to use as the source
                        (a Bio::EnsEMBL::Registry alias, defined in config file)

   --to string          The target species.
                        (a Bio::EnsEMBL::Registry alias, defined in config file)
                        More than one target species can be specified by using
                        several --to arguments or a comma-separated list, e.g.
                        -from human -to dog -to opossum or --to dog,opossum

   --release            The current Ensembl release. Needed for projection_info
                        database.

  [--compara string]    A Compara database
                        (a Bio::EnsEMBL::Registry alias, defined in config file)
                        If not specified, the first compara database defined in
                        the registry file is used.

  [--names]             Project display names and descriptions. Note only HGNC
                        display names are projected by default (see --all_sources).

  [--go_terms]          Project GO terms.

  [--delete_names]      Delete projected display xrefs & gene names.

  [--delete_go_terms]   Delete projected GO terms.

  [--print]             Print details of projection only, don't store in database

  [--method]            Type of homologs (default: TREE_HOMOLOGIES)

  [--all_sources]       Use all sources for name projection (default is HGNC only).

  [--nobackup]          Skip dumping of table backups

  [--full_stats]        Print full statistics, i.e.:
                         - number of terms per evidence type for projected GO terms
                         - number of display names per source (if using --all_sources)

  [--no_database]       Don't write to the projection_info database.

  [--quiet]             Don't print percentage progress information to STDOUT.

  [--one_to_many]       Also project one-to-many orthologs; multiple orthologs in
                        target are named "1 of 3", etc. Currently only affects
                        display xrefs, not GO terms.

  [--go_check]          Check if GO term is already assigned, and don't project if it
                        is. Off by default.

  [--help]              This text.

  Note that projected names or GO terms can be deleted from a database without doing any subsequent 
  projection by specifying only the -to, -delete_only and -delete_go_terms or -delete_names options.

  e.g

  perl project_display_xrefs.pl --conf compara_only.ini --host HOST -user USER -pass PASS -version 47 -names -delete_names -from human -to dog -nobackup -no_database

EOF

}

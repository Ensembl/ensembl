use strict;

# Sets display_xref_ids for novel genes in the "to" database based
# on their orthologs in the "from" database. Can also project GO xrefs.
# Orthology relationships are read from a Compara database.

use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;

my $method_link_type = "ENSEMBL_ORTHOLOGUES";

my ($conf, $compara, $from_species, @to_multi, $print, $names, $go_terms, $delete_names, $delete_go_terms);

GetOptions('conf=s'          => \$conf,
	   'compara=s'       => \$compara,
	   'from=s'          => \$from_species,
	   'to=s'            => \@to_multi,
	   'method=s'        => \$method_link_type,
	   'names'           => \$names,
	   'go_terms'        => \$go_terms,
	   'print'           => \$print,
	   'delete_names'    => \$delete_names,
	   'delete_go_terms' => \$delete_go_terms,
	   'help'            => sub { usage(); exit(0); });

@to_multi = split(/,/,join(',',@to_multi));

if (!$conf|| !$compara || !$from_species || !@to_multi) {

  usage();
  exit(1);

}

if (!$go_terms && !$names) {

  print "One or both of --names or --go_terms must be specified\n";
  print "Use --help for more detailed usage informaion\n";
  exit(1);

}

# Take values from ENSEMBL_REGISTRY environment variable or from ~/.ensembl_init
# if no reg_conf file is given.
Bio::EnsEMBL::Registry->no_version_check(1);
Bio::EnsEMBL::Registry->load_all($conf);

my $mlssa = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'MethodLinkSpeciesSet');
my $ha    = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'Homology');
my $ma    = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'Member');

my $from_ga = Bio::EnsEMBL::Registry->get_adaptor($from_species, 'core', 'Gene');

foreach my $to_species (@to_multi) {

  my $to_ga   = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'Gene');
  my $to_dbea = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'DBEntry');

  delete_names() if ($delete_names);
  delete_go_terms() if ($delete_go_terms);

  my $mlss = $mlssa->fetch_by_method_link_type_registry_aliases($method_link_type, [$from_species, $to_species]);

  my $str = "gene display_xrefs and descriptions" if ($names);
  $str .= " and " if ($names && $go_terms);
  $str .= "GO terms" if ($go_terms);

  print "Projecting $str from $from_species to $to_species\n";

  # build hash of external db name -> ensembl object type mappings
  my %db_to_type = build_db_to_type($to_ga);

  print "$to_species, before projection: \n" . &get_stats($to_ga);

  # Get all genes, find homologies, set xrefs

  foreach my $gene_id (@{$from_ga->list_dbIDs}) {

    my $gene = $from_ga->fetch_by_dbID($gene_id);
    # next unless ($gene->biotype eq "protein_coding");

    my $member = $ma->fetch_by_source_stable_id("ENSEMBLGENE",$gene->stable_id);
    next unless (defined $member);

    my $homologies = $ha->fetch_all_by_Member_MethodLinkSpeciesSet($member, $mlss);

    project_homologies($homologies, $to_ga, $to_dbea, $names, $go_terms, $ma, %db_to_type);

  }

  print "$to_species, after projection: \n" . &get_stats($to_ga);

}

# ----------------------------------------------------------------------

# Project homologies from first to subsequent species
sub project_homologies() {

  my ($homologies, $to_ga, $to_dbea, $names, $go_terms, $ma, %db_to_type) = @_;

  foreach my $homology (@{$homologies}) {

    my @mas = @{$homology->get_all_Member_Attribute};
    my ($from_member, $from_attribute) = @{$mas[0]};
    my ($to_member, $to_attribute) = @{$mas[1]};

    # ----------------------------------------
    # Display names and descriptions

    project_display_names($to_ga, $to_dbea, $ma, $from_member, $to_member, %db_to_type) if ($names);

    # ----------------------------------------
    # GO terms

    project_go_terms($to_ga, $to_dbea, $ma, $from_attribute, $to_attribute) if ($go_terms);

    # ----------------------------------------

  }

}

# --------------------------------------------------------------------------------

sub project_display_names {

  my ($to_ga, $to_dbea, $ma, $from_member, $to_member, %db_to_type) = @_;

  my $to_gene = $to_ga->fetch_by_stable_id($to_member->stable_id());

  # if no display name set, do the projection
  if (!$to_gene->external_name()) {

    my $from_gene = $from_ga->fetch_by_stable_id($from_member->stable_id());
    my $dbEntry = $from_gene->display_xref();

    # TODO only do this for certain types of DBEntry?

    if ($dbEntry) {

      # Modify the dbEntry to indicate it's not from this species
      my $txt = " [from $from_species gene " . $from_gene->stable_id() . "]";
      $dbEntry->display_id($dbEntry->display_id() . $txt);

      # Add the xref to the "to" gene, or transcript or translation depending on what the
      # other xrefs from this dbname as assigned to (see build_db_to_type)
      # Note that if type is not found, it means that we're dealing with a db that has no
      # xrefs in the target database, e.g. MarkerSymbol in mouse -> rat
      # In this case just assign to transcripts

      my @to_transcripts = @{$to_gene->get_all_Transcripts};
      my $to_transcript = $to_transcripts[0];

      my $type = $db_to_type{$dbEntry->dbname()};

      if ($type eq "Gene") {

	$to_gene->add_DBEntry($dbEntry);
	$to_dbea->store($dbEntry, $to_gene, 'Gene') if (!$print);

      } elsif ($type eq "Transcript" || !$type) {
	
	$to_transcript->add_DBEntry($dbEntry);
	$to_dbea->store($dbEntry, $to_transcript, 'Transcript') if (!$print);

      } elsif ($type eq "Translation") {

	my $to_translation = $to_transcript->translation();
	return if (!$to_translation);
	$to_translation->add_DBEntry($dbEntry);
	$to_dbea->store($dbEntry, $to_translation, 'Translation') if (!$print);

      } else {

	warn("Can't deal with xrefs assigned to $type (dbname=" . $dbEntry->dbname . ")\n");
	return;

      }

      # Set gene status to "KNOWN", modify the description and update display_xref

      $to_gene->status("KNOWN");
      $to_gene->description($from_gene->description() . $txt) if ($from_gene->description());
      $to_gene->display_xref($dbEntry);
	
      print $to_gene->stable_id() . " --> " . $dbEntry->display_id() . "\n" if ($print);

      # update the gene so that the display_xref_id is set
      $to_ga->update($to_gene) if (!$print);

    }

  }

}

# --------------------------------------------------------------------------------

sub project_go_terms {

  my ($to_ga, $to_dbea, $ma, $from_attribute, $to_attribute) = @_;

  # GO xrefs are linked to translations, not genes
  my $from_translation = $ma->fetch_by_dbID($from_attribute->peptide_member_id())->get_Translation();
  my $to_translation   = $ma->fetch_by_dbID($to_attribute->peptide_member_id())->get_Translation();

  my $to_go_xrefs = $to_translation->get_all_DBEntries();

  DBENTRY: foreach my $dbEntry (@{$from_translation->get_all_DBEntries()}) {

    next if ($dbEntry->dbname() ne "GO" || !$dbEntry);

    # only project GO terms with non-IEA evidence codes
    foreach my $et (@{$dbEntry->get_all_linkage_types}){
      next DBENTRY if ($et eq "IEA");
    }

    # check that each from GO term isn't already projected
    next if go_xref_exists($dbEntry, $to_go_xrefs);

    # add linkage_type for projection to IEA (in the absence of a specific one for projections)
    #$dbEntry->flush_linkage_types();
    $dbEntry->add_linkage_type("IEA");

    my $txt = " [from $from_species translation " . $from_translation->stable_id() . "]";

    $dbEntry->display_id($dbEntry->display_id() . $txt);

    $to_translation->add_DBEntry($dbEntry);

    print $to_translation->stable_id() . " --> " . $dbEntry->display_id() . "\n" if ($print);

    $to_dbea->store($dbEntry, $to_translation, 'Translation') if (!$print);
    #print "stored " . $to_translation->stable_id() . " " . $dbEntry->display_id() . "\n";

  }

}

# --------------------------------------------------------------------------------

sub go_xref_exists {

  my ($dbEntry, $to_go_xrefs) = @_;

  foreach my $xref (@{$to_go_xrefs}) {

    if ($xref->dbname eq $dbEntry->dbname() &&
	$xref->primary_id eq $dbEntry->primary_id()) {
      return 1;
    }

  }

  return 0;

}

# ----------------------------------------------------------------------

sub get_stats {

  my ($to_ga) = @_;

  my $total_genes = count_rows($to_ga, "SELECT COUNT(*) FROM gene g");

  my $str = "Total genes: $total_genes\t";

  my $count = count_rows($to_ga, "SELECT COUNT(*) FROM gene g, xref x WHERE g.display_xref_id=x.xref_id AND g.display_xref_id IS NOT NULL AND x.display_label NOT LIKE '%[from%'");
  $str .= sprintf("Gene names: unprojected %d (%3.1f\%)" , $count, (100 * $count / $total_genes));

  $count = count_rows($to_ga, "SELECT COUNT(*) FROM gene g, xref x WHERE g.display_xref_id=x.xref_id AND x.display_label LIKE '%[from%'");
  $str .= sprintf(" projected %d (%3.1f\%)" , $count, (100 * $count / $total_genes));

  $count = count_rows($to_ga, "SELECT COUNT(*) FROM gene g WHERE display_xref_id IS NOT NULL");
  $str .= sprintf(" total genes with names %d (%3.1f\%)" , $count, (100 * $count / $total_genes));

  $str .= "\nGO xrefs: total ";
  $str .= &count_rows($to_ga, "SELECT COUNT(*) FROM xref x, external_db e WHERE e.external_db_id=x.external_db_id AND e.db_name='GO'");

  $str .= " projected ";
  $str .= &count_rows($to_ga, "SELECT COUNT(*) FROM xref x, external_db e WHERE e.external_db_id=x.external_db_id AND e.db_name='GO' AND x.display_label LIKE '%[from%'");

  $str .= "\n";

  return $str;

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

  print "Setting gene display_xrefs that were projected to NULL\n";
  my $sth = $to_ga->dbc()->prepare("UPDATE gene, xref SET gene.display_xref_id = null WHERE gene.display_xref_id=xref.xref_id AND xref.display_label LIKE '%[from%'");
  $sth->execute();

  print "Deleting projected xrefs and object_xrefs\n";
  $sth = $to_ga->dbc()->prepare("DELETE x, ox FROM xref x, object_xref ox WHERE x.xref_id=ox.xref_id AND x.display_label LIKE '%[from%'");
  $sth->execute();

}

# ----------------------------------------------------------------------

sub delete_go_terms {

  my ($to_ga) = @_;

  print "Deleting projected GO terms\n";

  my $sth = $to_ga->dbc()->prepare("DELETE x, ox FROM xref x, external_db e, object_xref ox WHERE x.xref_id=ox.xref_id AND x.external_db_id=e.external_db_id AND e.db_name='GO' AND x.display_label like '%[from%';");
  $sth->execute();

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

   --compara string     A Compara database
                        (a Bio::EnsEMBL::Registry alias, defined in config file)

   --from string        The species to use as the source
                        (a Bio::EnsEMBL::Registry alias, defined in config file)

   --to string          The target species.
                        (a Bio::EnsEMBL::Registry alias, defined in config file)
                        More than one target species can be specified by using
                        several --to arguments or a comma-separated list, e.g.
                        -from human -to dog -to opossum or --to dog,opossum

  [--names]             Project display names and descriptions.

  [--go_terms]          Project GO terms.

  [--delete_names]      Delete projected display xrefs & gene names.

  [--delete_go_terms]   Delete projected GO terms.

  [--print]             Print details of projection only, don't store in database

  [--method]            Type of homologs (default: ENSEMBL_ORTHOLOGUES)

  [--help]              This text.

  e.g

  perl project_display_xrefs.pl -conf reg_conf.pl -compara compara35 -from human -to dog

EOF

}

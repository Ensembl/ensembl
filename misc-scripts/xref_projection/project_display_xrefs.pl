use strict;

# Sets display_xref_ids for novel genes in the "to" database based
# on their orthologs in the "from" database. Can also project GO xrefs.
# Orthology relationships are read from a Compara database.

use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;

my $method_link_type = "ENSEMBL_ORTHOLOGUES";

my ($conf, $compara, $from_species, @to_multi, $print, $store, $names, $go_terms);

GetOptions('conf=s'       => \$conf,
	   'compara=s'    => \$compara,
	   'from=s'       => \$from_species,
	   'to=s'         => \@to_multi,
	   'method=s'     => \$method_link_type,
	   'names'        => \$names,
	   'go_terms'     => \$go_terms,
	   'print'        => \$print,
           'store'        => \$store,
	   'help'         => sub { usage(); exit(0); });

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

my $changed_names = 0;
my $changed_go_xrefs = 0;

foreach my $to_species (@to_multi) {

  my $str = "gene display_xrefs and descriptions" if ($names);
  $str .= " and " if ($names && $go_terms);
  $str .= "GO terms" if ($go_terms);

  print "Projecting $str from $from_species to $to_species\n";

  my $to_ga   = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'Gene');
  my $to_dbea = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'DBEntry');

  my $mlss = $mlssa->fetch_by_method_link_type_registry_aliases($method_link_type, [$from_species, $to_species]);

  # Get all genes, find homologies, set xrefs

 
  foreach my $gene_id (@{$from_ga->list_dbIDs}) {

    my $gene = $from_ga->fetch_by_dbID($gene_id);
    # next unless ($gene->biotype eq "protein_coding");

    my $member = $ma->fetch_by_source_stable_id("ENSEMBLGENE",$gene->stable_id);
    next unless (defined $member);

    my $homologies = $ha->fetch_all_by_Member_MethodLinkSpeciesSet($member, $mlss);

    project_homologies($homologies, $to_ga, $to_dbea, $names, $go_terms, $ma);

  }

  print "Projected $changed_names gene names, out of a total of " . scalar(@{$from_ga->list_dbIDs}) . "\n" if ($names);
  print "Projected $changed_go_xrefs GO terms\n" if ($go_terms);

}

# ----------------------------------------------------------------------

# Project homologies from first to subsequent species
sub project_homologies() {

  my ($homologies, $to_ga, $to_dbea, $names, $go_terms, $ma) = @_;

  my $changed = 0;

  foreach my $homology (@{$homologies}) {

    my @mas = @{$homology->get_all_Member_Attribute};
    my ($from_member, $from_attribute) = @{$mas[0]};
    my ($to_member, $to_attribute) = @{$mas[1]};

    # ----------------------------------------
    # Display names and descriptions

    if ($names) {

      $changed_names += project_display_names($to_ga, $to_dbea, $ma, $from_member, $to_member);

    }

    # ----------------------------------------
    # GO terms

    if ($go_terms) {

      $changed_go_xrefs += project_go_terms($to_ga, $to_dbea, $ma, $from_attribute, $to_attribute);

    }


    # ----------------------------------------

  }

  return $changed;

}

# --------------------------------------------------------------------------------

sub project_display_names {

  my ($to_ga, $to_dbea, $ma, $from_member, $to_member) = @_;

  my $to_gene = $to_ga->fetch_by_stable_id($to_member->stable_id());

  my $changed = 0;

  # if no display name set, do the projection
  if (!$to_gene->external_name()) {

    my $from_gene = $from_ga->fetch_by_stable_id($from_member->stable_id());
    my $dbEntry = $from_gene->display_xref();

    # TODO only do this for certain types of DBEntry?

    if ($dbEntry) {

      # Modify the dbEntry to indicate it's not from this species
      my $txt = " [from $from_species gene " . $from_gene->stable_id() . "]";
      $dbEntry->display_id($dbEntry->display_id() . $txt);

      # Add the xref to the "to" gene, set its status to "KNOWN", modify the description and update display_xref
      $to_gene->add_DBEntry($dbEntry);
      $to_gene->status("KNOWN");
      $to_gene->description($from_gene->description() . $txt) if ($from_gene->description());
      $to_gene->display_xref($dbEntry);
	
      print $to_gene->stable_id() . " --> " . $dbEntry->display_id() . "\n" if ($print);

      # store - need to add the DBEntry and also update the gene display xref id
      if ($store) {
	$to_dbea->store($dbEntry, $to_gene, 'Gene');
	$to_ga->update($to_gene);
      }

      $changed++;

    }

  }

  return $changed;

}

# --------------------------------------------------------------------------------

sub project_go_terms {

  my ($to_ga, $to_dbea, $ma, $from_attribute, $to_attribute) = @_;

  my $changed = 0;

  # GO xrefs are linked to translations, not genes
  my $from_translation = $ma->fetch_by_dbID($from_attribute->peptide_member_id())->get_Translation();
  my $to_translation   = $ma->fetch_by_dbID($to_attribute->peptide_member_id())->get_Translation();

  my $to_go_xrefs = $to_translation->get_all_DBEntries();

  foreach my $dbEntry (@{$from_translation->get_all_DBEntries()}) {

    next if $dbEntry->dbname() ne "GO";

    # check that each from GO term isn't already projected, then project it
    next if go_xref_exists($dbEntry, $to_go_xrefs);

    # linkage types etc should be taken care of automatically
    $to_translation->add_DBEntry($dbEntry);

    my $txt = " [from $from_species translation " . $from_translation->stable_id() . "]";

    $dbEntry->display_id($dbEntry->display_id() . $txt);

    print $to_translation->stable_id() . " --> " . $dbEntry->display_id() . "\n" if ($print);

    $to_dbea->store($dbEntry, $to_translation, 'Translation') if ($store);

    $changed++;

  }

  return $changed;

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

sub usage {

  print << "EOF";

  Sets display_xref_ids and/or GO terms for novel genes in the "to" database
  based on their orthologs in the "from" database. Orthology relationships
  are read from a Compara database.

 perl project_display_xrefs.pl {options}

 Options ([..] indcicates optional):

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

   --names              Project display names and descriptions.

   --go_terms           Project GO terms.

  [--print]             Print details of projection

  [--store]             Upload projections to target databases

  [--method]            Type of homologs (default: ENSEMBL_ORTHOLOGUES)

  [--help]              This text.

  e.g

  perl project_display_xrefs.pl -conf reg_conf.pl -compara compara35 -from human -to dog

EOF

}

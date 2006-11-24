use strict;

# Sets display_xref_ids for novel genes in the "to" database based
# on their orthologs in the "from" database. Can also project GO xrefs.
# Orthology relationships are read from a Compara database.

use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;

my $method_link_type = "ENSEMBL_ORTHOLOGUES";

my ($conf, $compara, $from_species, @to_multi, $print, $names, $go_terms, $delete_names, $delete_go_terms, $no_backup, $full_stats, $descriptions, $release, $no_database);

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
	   'nobackup'        => \$no_backup,
	   'full_stats'      => \$full_stats,
           'descriptions'    => \$descriptions,
	   'release=i'       => \$release,
	   'no_database'     => \$no_database,
	   'help'            => sub { usage(); exit(0); });


if (!$conf) {

  print STDERR "Configuration file must be supplied via -conf argument\n";
  usage();
  exit(1);

} elsif (!$compara) {

  print STDERR "Compara database must be supplied via -compara argument\n";
  usage();
  exit(1);

} elsif (!$from_species) {

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

if (!$go_terms && !$names) {

  print STDERR "One or both of --names or --go_terms must be specified\n";
  print STDERR "Use --help for more detailed usage informaion\n";
  exit(1);

}

# only these evidence codes will be considered for GO term projection
my @evidence_codes = ( "IDA", "IEP", "IGI", "IMP", "IPI" );

@to_multi = split(/,/,join(',',@to_multi));

# Take values from ENSEMBL_REGISTRY environment variable or from ~/.ensembl_init
# if no reg_conf file is given.
Bio::EnsEMBL::Registry->no_version_check(1);
Bio::EnsEMBL::Registry->load_all($conf);

my $mlssa = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'MethodLinkSpeciesSet');
my $ha    = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'Homology');
my $ma    = Bio::EnsEMBL::Registry->get_adaptor($compara, 'compara', 'Member');

my $from_ga = Bio::EnsEMBL::Registry->get_adaptor($from_species, 'core', 'Gene');

my %projections_by_evidence_type;

foreach my $to_species (@to_multi) {

  my $to_ga   = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'Gene');
  my $to_dbea = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'DBEntry');

  write_to_projection_db($to_ga->dbc(), $release, $from_species, $from_ga->dbc(), $to_species) unless ($no_database);

  backup($to_ga) if (!$no_backup);

  delete_names($to_ga) if ($delete_names);
  delete_go_terms($to_ga) if ($delete_go_terms);

  my $mlss = $mlssa->fetch_by_method_link_type_registry_aliases($method_link_type, [$from_species, $to_species]);

  my $str = "gene display_xrefs and descriptions" if ($names);
  $str .= " and " if ($names && $go_terms);
  $str .= "GO terms" if ($go_terms);

  print "Projecting $str from $from_species to $to_species\n";

  # build hash of external db name -> ensembl object type mappings
  my %db_to_type = build_db_to_type($to_ga);

  print "$to_species, before projection: \n";
  print_stats($to_ga);

  # Get all genes, find homologies, set xrefs
  my @genes = @{$from_ga->list_dbIDs};
  my $i = 0;
  my $total_genes = scalar(@genes);
  foreach my $gene_id (@genes) {

    my $gene = $from_ga->fetch_by_dbID($gene_id);
    # next unless ($gene->biotype eq "protein_coding");

    print "$i of $total_genes source genes\n" if ($i % 1000 == 0);
    $i++;

    my $member = $ma->fetch_by_source_stable_id("ENSEMBLGENE",$gene->stable_id);
    next unless (defined $member);

    my $homologies = $ha->fetch_all_by_Member_MethodLinkSpeciesSet($member, $mlss);

    project_homologies($homologies, $to_ga, $to_dbea, $names, $go_terms, $ma, %db_to_type);

  }

  print "$to_species, after projection: \n";
  print_stats($to_ga);

  # print statistics if required
  if ($full_stats) {
    print "Projected terms by evidence code:\n";
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

}

# ----------------------------------------------------------------------

# Project homologies from first to subsequent species
sub project_homologies() {

  my ($homologies, $to_ga, $to_dbea, $names, $go_terms, $ma, %db_to_type) = @_;

  foreach my $homology (@{$homologies}) {

    next if ($homology->description() ne "ortholog_one2one" && $homology->description() ne "apparent_ortholog_one2one");

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
  my $from_gene = $from_ga->fetch_by_stable_id($from_member->stable_id());
  my $dbEntry = $from_gene->display_xref();
  my $to_source = $to_gene->display_xref()->dbname() if ($to_gene->display_xref());
  my $from_source = $from_gene->display_xref()->dbname() if ($from_gene->display_xref());

  my $from_latin_species = ucfirst(Bio::EnsEMBL::Registry->get_alias($from_species));

  # if no display name set, do the projection
  if (check_overwrite_display_xref($to_gene, $from_source, $to_source)) {

    if ($dbEntry) {

      # Modify the dbEntry to indicate it's not from this species - set info_type & info_text
      my $txt = "from $from_latin_species gene " . $from_gene->stable_id();

      # Add description to the gene if required
      # This should probably be in another column in the xref table but we just use the gene
      # description column for now.
      $to_gene->description($from_gene->description()) if ($descriptions && $from_gene->description());

      $dbEntry->info_type("PROJECTION");
      $dbEntry->info_text($txt);

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
	$to_dbea->store($dbEntry, $to_gene->dbID(), 'Gene') if (!$print);

      } elsif ($type eq "Transcript" || !$type) {
	
	$to_transcript->add_DBEntry($dbEntry);
	$to_dbea->store($dbEntry, $to_transcript->dbID(), 'Transcript') if (!$print);

      } elsif ($type eq "Translation") {

	my $to_translation = $to_transcript->translation();
	return if (!$to_translation);
	$to_translation->add_DBEntry($dbEntry);
	$to_dbea->store($dbEntry, $to_translation->dbID(), 'Translation') if (!$print);

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

    }

  }

}

# --------------------------------------------------------------------------------

sub project_go_terms {

  my ($to_ga, $to_dbea, $ma, $from_attribute, $to_attribute) = @_;

  # GO xrefs are linked to translations, not genes
  my $from_translation = $ma->fetch_by_dbID($from_attribute->peptide_member_id())->get_Translation();
  my $to_translation   = $ma->fetch_by_dbID($to_attribute->peptide_member_id())->get_Translation();

  my $from_latin_species = ucfirst(Bio::EnsEMBL::Registry->get_alias($from_species));

  my $to_go_xrefs = $to_translation->get_all_DBEntries();

  DBENTRY: foreach my $dbEntry (@{$from_translation->get_all_DBEntries()}) {

    next if ($dbEntry->dbname() ne "GO" || !$dbEntry);

    # only project GO terms with non-IEA evidence codes
    # also exclude ISS terms (manually projected based on orthologs)
    next if (ref($dbEntry) ne "Bio::EnsEMBL::GoXref");

    # TODO - this will skip whole xref if any evidence type is IEA
    # even if there are more than one evidence type for this GO term
    # Should be changed to just skip IEA one, not others
    foreach my $et (@{$dbEntry->get_all_linkage_types}){
      print "$et " . "\n";
      next DBENTRY if (!grep(/$et/, @evidence_codes));
    }

    # check that each from GO term isn't already projected
    next if go_xref_exists($dbEntry, $to_go_xrefs);

    # record statistics by evidence type
    foreach my $et (@{$dbEntry->get_all_linkage_types}){
      $projections_by_evidence_type{$et}++;
      #print $dbEntry->display_id() . " " . $et . " " . $projections_by_evidence_type{$et} . "\n";
    }

    # Change linkage_type for projection to IEA (in the absence of a specific one for projections)
    $dbEntry->flush_linkage_types();
    $dbEntry->add_linkage_type("IEA");

    my $txt = "from $from_latin_species translation " . $from_translation->stable_id();
    $dbEntry->info_type("PROJECTION");
    $dbEntry->info_text($txt);

    $to_translation->add_DBEntry($dbEntry);

    print $to_translation->stable_id() . " --> " . $dbEntry->display_id() . "\n" if ($print);

    $to_dbea->store($dbEntry, $to_translation->dbID(), 'Translation') if (!$print);
    #print "stored xref ID " . $dbEntry->dbID() ." " . $to_translation->stable_id() . " ". $to_translation->dbID() . " " . $dbEntry->display_id() . "\n";

  }

}

# --------------------------------------------------------------------------------

sub go_xref_exists {

  my ($dbEntry, $to_go_xrefs) = @_;

  foreach my $xref (@{$to_go_xrefs}) {

    next if (ref($dbEntry) ne "Bio::EnsEMBL::GoXref" || ref($xref) ne "Bio::EnsEMBL::GoXref");

    if ($xref->dbname() eq $dbEntry->dbname() &&
	$xref->primary_id() eq $dbEntry->primary_id() &&
	join("", @{$xref->get_all_linkage_types()}) eq join("", @{$dbEntry->get_all_linkage_types()})) {
      return 1;
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

    $count = count_rows($to_ga, "SELECT COUNT(*) FROM gene g, xref x WHERE g.display_xref_id=x.xref_id AND x.info_type='PROJECTION'");
    printf(" projected %d (%3.1f\%)" , $count, (100 * $count / $total_genes));

    $count = count_rows($to_ga, "SELECT COUNT(*) FROM gene g, xref x, external_db e WHERE g.display_xref_id=x.xref_id AND x.external_db_id=e.external_db_id AND e.db_name IN ('RefSeq_dna_predicted', 'RefSeq_peptide_predicted')");
    printf(" predicted %d (%3.1f\%)" , $count, (100 * $count / $total_genes));

    $count = count_rows($to_ga, "SELECT COUNT(*) FROM gene g WHERE display_xref_id IS NOT NULL");
    printf(" total genes with names %d (%3.1f\%)\n" , $count, (100 * $count / $total_genes));

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

  # do both old style (where display_label was modified) and new style (where info_type=PROJECTION)
  print "Setting gene display_xrefs that were projected to NULL\n";
  my $sth = $to_ga->dbc()->prepare("UPDATE gene, xref SET gene.display_xref_id = null WHERE gene.display_xref_id=xref.xref_id AND xref.info_type='PROJECTION'");
  $sth->execute();

  print "Deleting projected xrefs and object_xrefs\n";
  $sth = $to_ga->dbc()->prepare("DELETE x, ox FROM xref x, object_xref ox WHERE x.xref_id=ox.xref_id AND x.info_type='PROJECTION'");
  $sth->execute();

}

# ----------------------------------------------------------------------

sub delete_go_terms {

  my ($to_ga) = @_;

  print "Deleting projected GO terms\n";

  # do both old style (where display_label was modified) and new style (where info_type=PROJECTION)
  my $sth = $to_ga->dbc()->prepare("DELETE x, ox, gx FROM xref x, external_db e, object_xref ox, go_xref gx WHERE x.xref_id=ox.xref_id AND x.external_db_id=e.external_db_id AND ox.object_xref_id=gx.object_xref_id AND e.db_name='GO' AND x.info_type='PROJECTION'");
  $sth->execute();

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

  my ($to_gene, $from_dbname, $to_dbname) = @_;

  return 1 if (!$to_gene->external_name());

  if ($to_dbname eq "RefSeq_dna_predicted" || $to_dbname eq "RefSeq_peptide_predicted") {

    if (($from_species eq "human" && $from_dbname eq "HUGO") ||
	($from_species eq "mouse" && $from_dbname eq "MarkerSymbol")) {

      return 1;

    }

  }

  return 0;

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
    unless (system("mysql -h$host -P$port -u$user -p$pass -N -e 'select * from $table' $dbname > $dbname.$table.backup; gzip -9 -f $dbname.$table.backup") == 0) {
      print STDERR "Can't dump the original $table from $dbname for backup\n";
      exit 1;
    } else {
      print STDERR "Original $table table backed up in $dbname.$table.backup\n";
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

  my $sth = $db->prepare("INSERT INTO projections (release, timestamp, from_db, from_species_latin, from_species_common, to_db, to_species_latin, to_species_common) VALUES (?, NOW(), ?, ?, ?, ?, ?, ?)");

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

   --release            The current Ensembl release. Needed for projection_info
                        database.

  [--names]             Project display names and descriptions.

  [--go_terms]          Project GO terms.

  [--delete_names]      Delete projected display xrefs & gene names.

  [--delete_go_terms]   Delete projected GO terms.

  [--descriptions]      Project descriptions as well. Only works if -names is
                        specified. Descriptions appended to info_text field,
                        separated by |

  [--print]             Print details of projection only, don't store in database

  [--method]            Type of homologs (default: TREE_HOMOLOGIES)

  [--nobackup]          Skip dumping of table backups

  [--full_stats]        Print full statistics, e.g. number of terms per evidence type

  [--no_database]       Don't write to the projection_info database.

  [--help]              This text.

  e.g

  perl project_display_xrefs.pl -conf reg_conf.pl -compara compara35 -from human -to dog

EOF

}

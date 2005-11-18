use strict;

# Sets display_xref_ids for novel genes in the "to" database based
# on their orthologs in the "from" database. Orthology relationships 
# are read from a Compara database.

use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;

my $method_link_type = "ENSEMBL_ORTHOLOGUES";

my ($conf, $compara, $from_species, $to_species, $print, $store);

GetOptions('conf=s'       => \$conf,
	   'compara=s'    => \$compara,
	   'from=s'       => \$from_species,
	   'to=s'         => \$to_species,
	   'method=s'     => \$method_link_type,
	   'print'        => \$print,
           'store'        => \$store,
	   'help'         => sub { usage(); exit(0); });

if (!$compara || !$from_species || !$to_species) {

  usage();
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
my $to_ga   = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'Gene');
my $to_dbea = Bio::EnsEMBL::Registry->get_adaptor($to_species, 'core', 'DBEntry');

my $mlss = $mlssa->fetch_by_method_link_type_registry_aliases($method_link_type, [$from_species, $to_species]);

# Get all genes, find homologies, set xrefs

my $changed = 0;
my $unchanged = 0;

foreach my $gene_id (@{$from_ga->list_dbIDs}) {

  my $gene = $from_ga->fetch_by_dbID($gene_id);
  next unless ($gene->biotype eq "protein_coding");

  my $member = $ma->fetch_by_source_stable_id("ENSEMBLGENE",$gene->stable_id);
  next unless (defined $member);

  my $homologies = $ha->fetch_all_by_Member_MethodLinkSpeciesSet($member, $mlss);

  project_homologies($homologies);

}

print "Modified $changed genes, didn't modify $unchanged genes\n";

# ----------------------------------------------------------------------

# Project homologies from first to subsequent species
sub project_homologies() {

  my $homologies = shift;

 foreach my $homology (@{$homologies}) {

   my @mas = @{$homology->get_all_Member_Attribute};
   my ($from_member, $from_attribute) = @{$mas[0]};
   my ($to_member, $to_attribute) = @{$mas[1]};

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

       # Add the xref to the "to" gene, set its status to "KNOWN", modify the description and update display_xref
       $to_gene->add_DBEntry($dbEntry);
       $to_gene->status("KNOWN");
       $to_gene->description($from_gene->description() . $txt) if ($from_gene->description());
       $to_gene->display_xref($dbEntry);

       if ($print) {
	 print $to_gene->stable_id() . " --> " . $dbEntry->display_id() . "\n";
       }

       # store - need to add the DBEntry and also update the gene display xref id
       if ($store) {
	 $to_dbea->store($dbEntry, $to_gene, 'Gene');
	 $to_ga->update($to_gene);
       }

       $changed++;
     }

   } else {

     $unchanged++;

   }

 }

}


# ----------------------------------------------------------------------

sub usage {

  print << "EOF";

  Sets display_xref_ids for novel genes in the "to" database based on their
  orthologs in the "from" database. Orthology relationships are read from a
  Compara database.

 perl project_display_xrefs.pl {options}

 Options ([..] indcicates optional):

  [--conf filepath]     the Bio::EnsEMBL::Registry configuration file. If none
                        is given, the one set in ENSEMBL_REGISTRY will be used
                        if defined, if not ~/.ensembl_init will be used.

   --compara string     A Compara database
                        (a Bio::EnsEMBL::Registry alias, defined in config file)

   --from string        The species to use as the source
                        (a Bio::EnsEMBL::Registry alias, defined in config file)

   --to string          The target species
                        (a Bio::EnsEMBL::Registry alias, defined in config file)

  [--print]             Print details of projection

  [--store]             Upload projections to target databases

  [--method]            Type of homologs (default: ENSEMBL_ORTHOLOGUES)

  [--help]              This text.

  e.g

  perl project_display_xrefs.pl -conf reg_conf.pl -compara compara35 -from human -to dog

EOF

}

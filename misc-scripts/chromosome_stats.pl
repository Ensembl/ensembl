# Generate stats about genomes/chromosomes

use strict;
use warnings;

use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

# ----------------------------------------------------------------------
# Configuraton - change the values below to suit

# Core database configuration
my $core_host = "127.0.0.1";
my $core_port = 5001;
my $core_user = "ensro";
my $core_password = "";
my $core_db = "homo_sapiens_core_18_34";

# Compara database configuration (for homolog and family statistics)
my $compara_host = "127.0.0.1";
my $compara_port = 5000;
my $compara_user = "ensro";
my $compara_password = "";
my $compara_db = "ensembl_compara_18_1";

my $chr_name = "1"; # name of chromosome to perform calculations for

my $verbose = 1;    # if > 0, print debug output preceeded by a #

# ----------------------------------------------------------------------
# No user-serviceable parts inside

debug("Connecting to CORE database $core_db");
my $core_dbi = DBI->connect("dbi:mysql:host=$core_host;port=$core_port;database=$core_db", "$core_user", "$core_password", 
			    {'RaiseError' => 1}) || die "Can't connect to core DB";

debug("Connecting to COMPARA database $compara_db");
my $compara_dbi = DBI->connect("dbi:mysql:host=$compara_host;port=$compara_port;database=$compara_db", "$compara_user", "$compara_password", 
			       {'RaiseError' => 1}) || die "Can't connect to core DB";

# ----------------------------------------
# Get known, novel, pseudogene counts

debug("Getting gene counts for chromosome " . $chr_name);

my $core_dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(-user   => $core_user,
						  -dbname => $core_db,
						  -host   => $core_host,
						  -port   => $core_port,
						  -pass   => $core_password,
						  -driver => 'mysql' );

my $sa = $core_dba->get_SliceAdaptor();

my $slice = $sa->fetch_by_chr_name($chr_name);
my @genes = @{$slice->get_all_Genes()};

report("Total genes:\t" . scalar(@genes));

debug("Calculating known, novel and pseudogenes for chromosome " . $chr_name);

my $known_genes = 0;
my $novel_genes = 0;
my $pseudogenes = 0;

foreach my $gene (@genes) {
  #print $gene->stable_id . " " . $gene->is_known() . " " . $gene->type() . "\n";
  $known_genes++ if $gene->is_known  && $gene->type eq 'ensembl';
  $novel_genes++ if !$gene->is_known && $gene->type eq 'ensembl';
  $pseudogenes++ if $gene->type eq 'pseudogene';
}

report("Known genes:\t" . $known_genes);
report("Novel genes:\t" . $novel_genes);
report("Pseudogenes:\t" . $pseudogenes);

# ----------------------------------------
# Get numbers of genes with homologs in other species

debug("Calculating homologs for chromosome " . $chr_name);

# the SQL below will calculate the number of homologous genes between
# human and all other species to use a "source" species other than
# human, modify m1.genome_db_id = 1 to whatever the compara
# genome_db_id is for the required species.

my $homolog_sql = "select count( distinct m1.stable_id ) as num_genes, g.name from homology_member hm1, homology_member hm2, member m1, member m2, genome_db g where hm1.member_id = m1.member_id and m1.chr_name = '" . $chr_name . "' and m1.genome_db_id = 1 and hm2.member_id = m2.member_id and hm1.homology_id = hm2.homology_id and hm1.member_id != hm2.member_id and m2.genome_db_id=g.genome_db_id GROUP BY g.name";

my $sth = $compara_dbi->prepare($homolog_sql);
my ($num_genes, $species);
$sth->execute();
$sth->bind_columns(\$num_genes, \$species);
report("Number of homologous genes between human and other species");
while (my @row = $sth->fetchrow_array()) {
  report($species . "\t" . $num_genes);
}
$sth->finish();

# ----------------------------------------
# How many genes are in a family (of proteins) that has one or more members in other species


# --------------------------------------------------------------------------------

sub debug {

  my $str = shift;

  print "# " . $str . "\n" if $verbose;

}

sub report {

  my $str = shift;

  print $str . "\n";

}

# ----------------------------------------------------------------------

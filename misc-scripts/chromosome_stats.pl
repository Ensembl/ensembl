# Generate stats about genomes/chromosomes

use strict;
use warnings;

use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

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

my $chr_name = "1";   # name of chromosome to perform calculations for

my %family_species = (6 => "drosophila", # set this to the species you want for family comparisons;
		      7 => "C. elegans"); # see the genome_db table in the compara database

my $family_coverage_percentage = 40.0; # families who have more than this %age of genes on $chr_name will be printed

my $exon_density_step_size = 100000; # step size for exon density table
my $exon_density_file_name = "exon_density_chr_" . $chr_name;

my $verbose = 1;	 # if > 0, print debug output preceeded by a #

# set the following to 1 to have stats printed
my $do_gene_counts = 1;
my $do_homology = 1;
my $do_family = 1;
my $do_exon_density = 1;

# ----------------------------------------------------------------------
# No user-serviceable parts inside

my $core_dbi;
if ($do_gene_counts) {
  debug("Connecting to CORE database $core_db");
  $core_dbi = DBI->connect("dbi:mysql:host=$core_host;port=$core_port;database=$core_db", "$core_user", "$core_password", {'RaiseError' => 1}) || die "Can't connect to core DB";
}

my $compara_dbi;
if ($do_homology || $do_family) {
  debug("Connecting to COMPARA database $compara_db");
  $compara_dbi = DBI->connect("dbi:mysql:host=$compara_host;port=$compara_port;database=$compara_db", "$compara_user", "$compara_password", {'RaiseError' => 1}) || die "Can't connect to core DB";
}

my %peptide_gene = (); # will store mapping between peptide ID and gene ID
my $sth;
my $slice;
my $rr;

# ----------------------------------------
# Get known, novel, pseudogene counts

if ($do_gene_counts || $do_family || $do_exon_density) { # need to build peptide_gene for family too

  debug("Getting gene counts for chromosome " . $chr_name);

  my $core_dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(-user   => $core_user,
						    -dbname => $core_db,
						    -host   => $core_host,
						    -port   => $core_port,
						    -pass   => $core_password,
						    -driver => 'mysql' );

  my $sa = $core_dba->get_SliceAdaptor();

  $slice = $sa->fetch_by_chr_name($chr_name);
  my @genes = @{$slice->get_all_Genes()};

  $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();
  foreach my $gene (@genes) {
    my $gene_stable_id = $gene->stable_id;
    my @transcripts = @{$gene->get_all_Transcripts()};
    foreach my $transcript (@transcripts) {
      my $translation = $transcript->translation;
      next if !$translation; # check with Ewan if pseudogenes are to be included
      my $translation_stable_id = $translation->stable_id;
      $peptide_gene{$translation_stable_id} = $gene_stable_id;

      if ($do_exon_density) {
	
	# Add each exon to the RangeRegistry
	my @exons = @{$transcript->get_all_Exons()};
	foreach my $exon (@exons) {
	  $rr->check_and_register($chr_name, $exon->start, $exon->end);
	}

      }
    }
  }

  report("Total genes:\t" . scalar(@genes));

  if ($do_gene_counts) {

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
  }

}

# ----------------------------------------
# Exon density
# Note all exons should have been checked and registered
# with the RangeRegistry for this to work

if ($do_exon_density) {

  debug("Caclulating exon density for steps of " . $exon_density_step_size . " bases on chromosome " . $chr_name);

  open (EXON_FILE, ">" . $exon_density_file_name);

  my %exon_densities = (); # will hold base pair index -> fractional exon density

  my $chunk_start = 1;
  my $chunk_end = $exon_density_step_size;

  while ($chunk_end < $slice->chr_end()) {

    my $non_exon_length = 0;

    my $result = $rr->check_and_register($chr_name, $chunk_start, $chunk_end);
    # result (if defined) is a listref of [start, end] pairs; in this case
    # these are the start and end of the "bits between exons"; so we can get
    # the number of bases in exons in this chunk by subtracting the total
    # from the length of the chunk
    if (defined $result) {
      my @bits = @{$result};
      foreach my $bit (@bits) {                # bit is an array reference
	my $len = ($bit->[1] - $bit->[0]) + 1;
	$non_exon_length += $len;
      }
    }

    # calculate exon_density for this chunk
    my $exon_density = ($exon_density_step_size - $non_exon_length) / $exon_density_step_size;
    print EXON_FILE $chunk_start . "->" . $chunk_end . "\t" . $exon_density;

    $chunk_start += $exon_density_step_size;
    $chunk_end += $exon_density_step_size;

  }

  # do final bit between $chunk_end and chr_end
  # TODO

  debug("Wrote exon density data to " . $exon_density_file_name);

}

# ----------------------------------------
# How many genes are in a family (of proteins) that has one or more members in other species

if ($do_family) {

  foreach my $target_genome_db_id (keys %family_species) {

    my $species_name = $family_species{$target_genome_db_id};

    debug("Getting family data for human / " . $species_name . ", human chromosome " . $chr_name);

    my %human_genes_in_other_species_with_family = ();

    my $family_sql = "select m1.stable_id, m1.genome_db_id, m1.chr_name, fm1.family_id from family_member fm1, member m1 where fm1.member_id = m1.member_id and m1.source_id = 1 order by fm1.family_id";

    $sth = $compara_dbi->prepare($family_sql);
    my ($peptide_id, $genome_db_id, $compara_chr_name, $family_id);
    $sth->execute();
    $sth->bind_columns(\$peptide_id, \$genome_db_id, \$compara_chr_name, \$family_id);
    my $current_family_id = -1;

    my %tmp_genes = ();
    my $other_species_seen = 0;
    my $genes_in_family = 0;

    while (my @row = $sth->fetchrow_array()) {

      # end of current family reached ?
      if ($family_id != $current_family_id) {

	if ($other_species_seen) {
	  # add to "global" list
	  foreach my $k (keys %tmp_genes) {
	    $human_genes_in_other_species_with_family{$k} = 1;
	  }
	}
	$other_species_seen = 0;

	# is this family of a reasonable (>4) size and are over 80%
	# of its genes on chr 1?
	if ($genes_in_family > 4) {
	  my $percentage = 100.0 * (scalar keys %tmp_genes) / $genes_in_family;
	  if ($percentage >= $family_coverage_percentage) {
	    report("Family with ID " . $current_family_id . " has " . $percentage ."% of its genes on chr " . $chr_name);
	  }
	}

	# reset counters
	$current_family_id = $family_id;
	%tmp_genes = ();
	$genes_in_family = 0;
      }

      $genes_in_family++;

      if ($compara_chr_name eq $chr_name) { # only interested in those from say chr 1

	# result is already ordered by family id
	# want to store human gene IDs corresponding to peptides,
	# *if* there is a drosophila (or whatever) peptide in the same family
	
	# store gene IDs of human genes
	if ($genome_db_id == 1) {
	  my $key = $peptide_gene{$peptide_id}; # Translate peptide IDs -> gene IDs
	  $tmp_genes{$key} = 1;
	}

      }

      # does this row correspond to the species we're interested in?
      if ($genome_db_id == $target_genome_db_id) {
	$other_species_seen = 1;
      }

    }

    $sth->finish();

    report("Number of human genes on chromosome " . $chr_name . " in families that also have members in " . $species_name . ": " . scalar keys %human_genes_in_other_species_with_family);

  }
}

# ----------------------------------------
# Get numbers of genes with homologs in other species

if ($do_homology) {

  debug("Calculating homologs for chromosome " . $chr_name);

  # the SQL below will calculate the number of homologous genes between
  # human and all other species to use a "source" species other than
  # human, modify m1.genome_db_id = 1 to whatever the compara
  # genome_db_id is for the required species.

  my $homolog_sql = "select count( distinct m1.stable_id ) as num_genes, g.name from homology_member hm1, homology_member hm2, member m1, member m2, genome_db g where hm1.member_id = m1.member_id and m1.chr_name = '" . $chr_name . "' and m1.genome_db_id = 1 and hm2.member_id = m2.member_id and hm1.homology_id = hm2.homology_id and hm1.member_id != hm2.member_id and m2.genome_db_id=g.genome_db_id GROUP BY g.name";

  $sth = $compara_dbi->prepare($homolog_sql);
  my ($num_genes, $species);
  $sth->execute();
  $sth->bind_columns(\$num_genes, \$species);
  report("Number of homologous genes between human and other species");
  while (my @row = $sth->fetchrow_array()) {
    report($species . "\t" . $num_genes);
  }
  $sth->finish();
}


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

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;



my $host = 'ecs1g';
my $user = '';
my $pass = '';
my $dbname = 'homo_sapiens_core_20_34';
my $port = '3306';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);


my $slice_adaptor = $db->get_SliceAdaptor();

my @chromosomes = @{$slice_adaptor->fetch_all('chromosome')};

foreach my $chromosome (@chromosomes) {
  my $num_known_genes  = 0;
  my $num_genes        = 0;
  my $num_pseudo_genes = 0;

  print "Processing chromosome ", $chromosome->seq_region_name(), "\n";

  my @genes = @{$chromosome->get_all_Genes()};

  foreach my $gene (@genes) {
    if($gene->type() eq 'pseudogene') {
      $num_pseudo_genes++;
    } else {
      $num_genes++;
      if($gene->is_known()) {
        $num_known_genes++;
      }
    }
  }

  $slice_adaptor->set_seq_region_attrib($chromosome,
                                        'GeneCount', $num_genes);
  $slice_adaptor->set_seq_region_attrib($chromosome,
                                        'KnownGeneCount', $num_known_genes);

  $slice_adaptor->set_seq_region_attrib($chromosome,
                                        'PseudoGeneCount', $num_pseudo_genes);
}


print_chromo_stats(\@chromosomes);

sub print_chromo_stats {
  my $chromosomes = shift;

  foreach my $chr (@$chromosomes) {
    print "\nchromosome: ",$chr->seq_region_name(),"\n";
    my ($num_genes)  = $chr->get_attribute('GeneCount');
    my ($num_known)  = $chr->get_attribute('KnownGeneCount');
    my ($num_pseudo) = $chr->get_attribute('PseudoGeneCount');
    print "NumGenes: $num_genes\n";
    print "NumKnown: $num_known\n";
    print "NumPseudo: $num_pseudo\n";
  }
}


1;



use strict;

BEGIN {
  my $ensembl_root = '/ensweb/wwwdev/server'; 
  unshift @INC, "$ensembl_root/modules";
  unshift @INC, "$ensembl_root/ensembl-draw/modules";
  unshift @INC, "$ensembl_root/ensembl/modules";
  unshift @INC, "$ensembl_root/ensembl-external/modules";
  unshift @INC, "$ensembl_root/ensembl-map/modules";
  unshift @INC, "$ensembl_root/ensembl-compara/modules";
  unshift @INC, "$ensembl_root/conf";
  unshift @INC, "$ensembl_root/bioperl-live";
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $host   = 'localhost';
my $user   = 'ensadmin';
my $pass   = 'ensembl';
my $dbname = 'mcvicker_homo_sapiens_core_20_34b';
my $port = '3306';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);

my $slice_adaptor = $db->get_SliceAdaptor();
my $attrib_adaptor = $db->get_AttributeAdaptor();

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

  my @attribs;

  push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Gene Count',
     -CODE => 'GeneCount',
     -VALUE => $num_genes,
     -DESCRIPTION => 'Total Number of Genes');

  push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Known Gene Count',
     -CODE => 'KnownGeneCount',
     -VALUE => $num_known_genes,
     -DESCRIPTION => 'Total Number of Known Genes');

  push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'PseudoGene Count',
     -CODE => 'PseudoGeneCount',
     -VALUE => $num_pseudo_genes,
     -DESCRIPTION => 'Total Number of PseudoGenes');

  $attrib_adaptor->store_on_Slice($chromosome, \@attribs);
}


print_chromo_stats(\@chromosomes);

sub print_chromo_stats {
  my $chromosomes = shift;

  foreach my $chr (@$chromosomes) {
    print "\nchromosome: ",$chr->seq_region_name(),"\n";
    foreach my $attrib (@{$chr->get_all_Attributes()}) {
      print "  ", $attrib->name(), ": ", $attrib->value(), "\n";
    }
  }
}


1;



use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname  );


GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname
	  );
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);

#
# Only run on database with genes
#
my $sth = $db->prepare( "select count(*) from gene" );
$sth->execute();

my ( $gene_count )  = $sth->fetchrow_array();

if( ! $gene_count ) {
  print STDERR "No gene density for $dbname.\n";
  exit();
}

#
# and seq_regions
#
$sth = $db->prepare( "select count(*) from seq_region" );
$sth->execute();
my ( $seq_region_count ) = $sth->fetchrow_array();
if( ! $seq_region_count ) {
  print STDERR "No seq_regions for $dbname.\n";
  exit();
}

my $slice_adaptor = $db->get_SliceAdaptor();
my $attrib_adaptor = $db->get_AttributeAdaptor();

my $top_slices = $slice_adaptor->fetch_all( "toplevel" );


foreach my $slice (@$top_slices) {
  my $num_known_genes  = 0;
  my $num_genes        = 0;
  my $num_pseudo_genes = 0;

  print "Processing seq_region ", $slice->seq_region_name(), "\n";

  my @genes = @{$slice->get_all_Genes()};

  foreach my $gene (@genes) {
    if(uc($gene->type()) eq uc('pseudogene')) {
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

  $attrib_adaptor->store_on_Slice($slice, \@attribs);
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



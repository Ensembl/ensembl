use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Lite::DBAdaptor;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname, $genestats, $snpstats  );


GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname,
      "genestats", \$genestats,
      "snpstats", \$snpstats
	  );
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);


# do both genestats and snpstats by default
$genestats = $snpstats = 1 if(!$genestats && !$snpstats);


#
# Only run on database with genes
#

my $genes_present;

if($genestats) {
  my $sth = $db->prepare( "select count(*) from gene" );
  $sth->execute();

  my ( $gene_count )  = $sth->fetchrow_array();

  $genes_present = ($gene_count) ? 1 : 0;
} else {
  $genes_present = 0;
}

#
# and seq_regions
#
my $sth = $db->prepare( "select count(*) from seq_region" );
$sth->execute();
my ( $seq_region_count ) = $sth->fetchrow_array();
if( ! $seq_region_count ) {
  print STDERR "No seq_regions for $dbname.\n";
  exit();
}

my $snps_present = $snpstats && lite_attach( $db );


my $slice_adaptor = $db->get_SliceAdaptor();
my $attrib_adaptor = $db->get_AttributeAdaptor();

my $top_slices = $slice_adaptor->fetch_all( "toplevel" );


foreach my $slice (@$top_slices) {
  print STDERR "Processing seq_region ", $slice->seq_region_name(), "\n";

  my @attribs;

  if($genes_present) {
    my $num_known_genes  = 0;
    my $num_genes        = 0;
    my $num_pseudo_genes = 0;

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
  }

  if( $snps_present ) {
    my $snps = $slice->get_all_SNPs();
    push @attribs, Bio::EnsEMBL::Attribute->new
      (-NAME => 'SNP Count',
       -CODE => 'SNPCount',
       -VALUE => scalar( @$snps ),
       -DESCRIPTION => 'Total Number of SNPs');
  }

  $attrib_adaptor->store_on_Slice($slice, \@attribs);
#  print_chromo_stats([$slice]);
}



sub print_chromo_stats {
  my $chromosomes = shift;

  foreach my $chr (@$chromosomes) {
    print "\nchromosome: ",$chr->seq_region_name(),"\n";
    foreach my $attrib (@{$chr->get_all_Attributes()}) {
      print "  ", $attrib->name(), ": ", $attrib->value(), "\n";
    }
  }
}


#
# tries to attach lite.
#

sub lite_attach {
  my $db = shift;

  my $core_db_name;
  $core_db_name = $db->dbname();
  if( $core_db_name !~ /_core_/ ) {
    return 0;
  }
  #
  # get a lost of all databases on that server
  #
  my $sth = $db->prepare( "show databases" );
  $sth->execute();
  my $all_db_names = $sth->fetchall_arrayref();
  my %all_db_names = map {( $_->[0] , 1)} @$all_db_names;
  my $snp_db_name = $core_db_name;
  $snp_db_name =~ s/_core_/_lite_/;
  if( ! exists $all_db_names{ $snp_db_name } ) {
    return 0;
  }

  my $snp_db = Bio::EnsEMBL::Lite::DBAdaptor->new
    ( -host => $db->host(),
      -user => $db->username(),
      -pass => $db->password(),
      -port => $db->port(),
      -dbname => $snp_db_name );
  $db->add_db_adaptor( "lite", $snp_db );
  return 1;
}


1;



#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config pass_through);

# (make sure api version is correct
# Usage:
# perl alt_alleles.pl -pass XXXX > & human_release_63_alt_alleles
#
#
# long way
# perl alt_alleles.pl -dbname homo_sapiens_core_63_37 -host ens-staging1 -pass XXXX > & human_release_63_alt_alleles
#

my ($host, $pass, $port, $dbname, $user);

GetOptions(
    'user=s'  => \$user,
    'pass=s'  => \$pass,
    'host=s'  => \$host,
    'port=i'  => \$port,
    'dbname=s'       => \$dbname);



my $api_version = Bio::EnsEMBL::ApiVersion->software_version();

if(!defined($dbname)){
  $dbname = "homo_sapiens_core_".$api_version."_37";
}

# Connect to the core database 

my $core_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host => $host||'ens-staging1',
					  -port => $port || '3306',
                                          -user => $user||'ensadmin',
					  -pass => $pass,
                                          -dbname => $dbname);

my $sa = $core_dba->get_adaptor("slice");
my $ga = $core_dba->get_adaptor("gene");

my @top_slices = @{$sa->fetch_all("toplevel",undef,1)}; #returns all slices including non-reference ("toplevel",undef,1)
print "Looping over chromosomes...\n";

my %alt_alleles;
my %reference_genes;

my %logic_names = ("havana" => 1,
		   "ensembl_havana_gene" => 1,
		   "ensembl_havana_lincrna" => 1,);

my %skip_sources = ("Clone_based_vega_gene" => 1,
		    "Clone_based_ensembl_gene" => 1);

foreach my $slice (@top_slices) {
  my $slice_name = $slice->seq_region_name;
 GENE:
  foreach my $gene (@{$ga->fetch_all_by_Slice($slice)}) {
     
    next GENE if (exists($skip_sources{$gene->external_db}) );

    next GENE unless (exists($logic_names{$gene->analysis->logic_name}) );

    my $gene_id = $gene->dbID;

    #get name and synonyms
    my ($syns,$display_id);
    eval { $display_id = $gene->display_xref->display_id }; 
    if ($@) {
      print "Can't get name for gene $gene_id on slice $slice_name\n";
    } else {
	push @{$alt_alleles{$display_id}}, $gene_id ;
	if ($slice->is_reference) {
	    $reference_genes{$gene_id} = 1;
	}

    }

  }

}   

#
# Delete the old data
#

my $sth = $core_dba->dbc->prepare("delete from alt_allele");
$sth->execute;


#
# Check that all alt_alleles are valid.
# NOTE: if an alt allele has only one gene_id delete it from the hash 
#

foreach my $key (keys %alt_alleles){
  my @arr = @{$alt_alleles{$key}};
  if(scalar(@arr) <= 1){
      delete $alt_alleles{$key};
  }
}

#
# Then store the alt_allele and gene_id in the alt_allele table.
#

my $sql = "insert into alt_allele (alt_allele_id, gene_id, is_ref) values(?, ?, ?)";
my $insert_sth = $core_dba->dbc->prepare($sql);

my $count=0;
my $alt_id_count = 0;



foreach my $key (keys %alt_alleles){
  my @arr = @{$alt_alleles{$key}};
  if(scalar(@arr) > 1){
    my $sanity_check =0; # should be 1 reference gene only if 0 or more than 1 we may have a problem
    $alt_id_count++;
    foreach my $gene_id (@arr){
      my $ref = 0;
      if(exists($reference_genes{$gene_id})){
	$ref = 1;
	$sanity_check++;
      }
      $insert_sth->execute($alt_id_count, $gene_id, $ref);
      $count++;
    }
    if($sanity_check !=1){
      print STDERR "Problem we have $sanity_check genes on the reference sequence for alt_id $alt_id_count gene display_id $key\n";
    }
   
  }

}

print "Added $alt_id_count alt_allele ids for $count genes\n";

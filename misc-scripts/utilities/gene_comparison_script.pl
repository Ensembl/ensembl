#!/usr/local/bin/perl


=head1 NAME

 gene_comparison_script.pl

=head1 DESCRIPTION

  This script is an example of use of the class Bio::EnsEMBL::Utils::GeneComparison
  retrieving data from two databases.
  In particular, this script gets some statistics on the whole and only-coding exon overlaps
  for two sets of genes. The default is to compare the ensembl build genes from homo_sapiens_core_110
  with human annotated genes on chr20.

  Notice that GeneComparison can be used with any sets of genes, not necessarily obtained in the way
  coded in this script. For more information, see the documentation of this class.

  WARNING: the default databases are in the old schema, so use branch code

=head1 SYNOPSIS

 We can take a piece of chromosome, e.g.

   gene_comparison_script.pl -chr chr20 -chrstart 3000000 -chrend 4000000

 we can specify the gene types we want to compare, e.g.
 
   gene_comparison_script.pl -chr chr20 -type1 HUMACE-Novel_CDS -type2 ensembl 
  
 also we can include more than two types, e.g.

   gene_comparison_script.pl -chr chr20 -type1 HUMACE-Novel_CDS -type1 HUMACE-Known -type2 ensembl 

 Other things that can be especified are: databases for each gene set to be compared ('dbname1' and 'dbname2')
 or one single database for both ('dbname'), golden-paths for each database ('path1' and 'path2') or a single one
 ('path'), the user ('user') and either one common host ('host') or one host for each 
 database ('host1' and 'host2'). See the GetOptions function below for more details												   
=cut 

use strict;  
use diagnostics;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::GeneComparison;
use Getopt::Long;


my $host1  = 'ecs1b';
my $host2  = 'ensrv3';
my $user   = 'ensro';
my $host;
my $dbname;
my $path;

# WARNING: these databases are in the old schema, so use branch code

# the default behaviour is to compare ensembl built genes with human-annotated genes in chr20

my $dbname1 = 'chr20';
my $dbname2 = 'homo_sapiens_core_110';
my $path1  = 'Sanger_02';
my $path2  = 'UCSC';
my $type1  = ['HUMACE-Novel_CDS','HUMACE-Known'];
my $type2  = ['ensembl'];
my (@type1,@type2);

# if one only provides $chr of these variables, the whole chromosome is taken #

my ($chrstart,$chrend);
my $chr;
my (@opt_type1,@opt_type2);


#change it so that using -type does not add up on top of the definitions above but

# actually substitute them

&GetOptions( 'host:s'    => \$host,
	     'host1:s'   => \$host1,
	     'host2:s'   => \$host2,
	     'dbuser:s'  => \$user,
	     'dbname:s'  => \$dbname,
	     'dbname1:s' => \$dbname1,
	     'dbname2:s' => \$dbname2,
	     'path:s'    => \$path,
	     'path1:s'   => \$path1,
	     'path2:s'   => \$path2,
	     'type1:s@'  => \@opt_type1, 
	     'type2:s@'  => \@opt_type2,
	     'chr=s'     => \$chr,
	     'chrstart:n'=> \$chrstart,
	     'chrend:n'  => \$chrend,
	     
	   );     
            
if ($host)  {
  $host1=$host;
  $host2=$host;
}

if ($dbname){
  $dbname1=$dbname;
  $dbname2=$dbname;
}

if ($path)  {
  $path1=$path;
  $path2=$path;
}

if (@opt_type1){
  $type1 = \@opt_type1;
}

if (@opt_type2){
  $type2 = \@opt_type2;
}


# connect to the database 
my $db1= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $host1,
					    -user  => $user,
					    -dbname=> $dbname1);

print STDERR "Connected to database $dbname1\n";

my $db2= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $host2,
					     -user  => $user,
					     -dbname=> $dbname2);

print STDERR "Connected to database $dbname2\n";

# use different golden paths
$db1->static_golden_path_type($path1); 
$db2->static_golden_path_type($path2); 

my $sgp1 = $db1->get_StaticGoldenPathAdaptor;
my $sgp2 = $db2->get_StaticGoldenPathAdaptor;

# get a virtual contig with a piece-of/entire chromosome #
my ($vcontig1,$vcontig2);
if ($chrstart && $chrend){
  print STDERR "Fetching region $chr, $chrstart - $chrend\n";
  $vcontig1 = $sgp1->fetch_VirtualContig_by_chr_start_end($chr,$chrstart,$chrend);
  $vcontig2 = $sgp2->fetch_VirtualContig_by_chr_start_end($chr,$chrstart,$chrend);
}
else{
  print STDERR "Fetching $chr\n";
  $vcontig1 = $sgp1->fetch_VirtualContig_by_chr_name($chr);
  $vcontig2 = $sgp2->fetch_VirtualContig_by_chr_name($chr);
}

# get the genes of type @type1 and @type2 from $vcontig1 and $vcontig2, respectively #
# the default in type1 is ensembl genes and in type2 HUMACE-Novel_CDS and HUMACE-Known #

my (@genes1,@genes2);

foreach my $type ( @{ $type1 } ){
  print STDERR "Fetching genes of type $type\n";
  my @more_genes = $vcontig1->get_Genes_by_Type($type);
  push ( @genes1, @more_genes ); 
  print STDERR scalar(@more_genes)." genes found\n";
}

foreach my $type ( @{ $type2 } ){
  print STDERR "Fetching genes of type $type\n";
  my @more_genes = $vcontig2->get_Genes_by_Type($type);
  push ( @genes2, @more_genes ); 
  print STDERR scalar(@more_genes)." genes found\n";
}

# get a GeneComparison object 
my $gene_comparison = Bio::EnsEMBL::Utils::GeneComparison->new(\@genes1, \@genes2);
# as convention, we put first the annotated (or benchmark) genes and second the predicted genes
# and the comparison methods refer to the second list with respect to the first one

## As an example, we get the number of exons per percentage overlap using coding exons only
#my %coding_statistics =  $gene_comparison->get_Coding_Exon_Statistics;

## You could also do it for all exons, coding and non-conding
#
# my %statistics = $gene_comparison->get_Exon_Statistics;
#

# The hashes hold the number of occurences as values and integer percentage overlap as keys
# these methods also print out the start and end coding exon overlaps

# We can produce a histogram

#my @values = values ( %coding_statistics );
#@values = sort {$b <=> $a} @values;

#print "Percentage overlap : Number of overlapping coding exons\n";
#for (my $i=1; $i<= 100; $i++){
#  if ( $coding_statistics{$i} ){
#    print $i." :\t".$coding_statistics{$i}."\t".print_row ($coding_statistics{$i})."\n";
#  }
#  else{
#    print $i." :\n";
#  }
#}

#sub print_row {
#  my $size = int( shift @_ );
#  $size = int( log( 1000*$size/($values[0]) ) ); # tweak this to re-scale it as you wish
#  my $row='';
#  for (my $i=0; $i<$size; $i++){
#    $row .='*';
#  }
#  return $row;
#}

#########################################################
#
#  Other examples of the potential use of GeneComparison
#
#########################################################

## cluster the genes we have passed to $gene_comparison

my @clusters    = $gene_comparison->cluster_Genes;

my @unclustered = $gene_comparison->unclustered_Genes;

# print the clusters 
print "Number of clusters: ".scalar( @clusters )."\n";

# my $count=1;
# foreach my $cluster (@clusters){
#   print "Cluster $count:\n";
#   print $cluster->to_String."\n";
#   $count++;
# }

#$count=1;
print "Unmatched genes: ".scalar( @unclustered )."\n";

#foreach my $cluster (@unclustered){
#print "Unclustered $count:\n";
#print $cluster->to_String."\n";
#$count++;
#}

# get the missing exons ( those in @genes1 which are missing in @genes2 )
$gene_comparison->find_missing_Exons(\@clusters);

# get the overpredicted exons ( those in @genes2 which are missing in @genes2 )
$gene_comparison->find_overpredicted_Exons(\@clusters);

## get the list of unmatched genes
# 
# my ($unmatched1,$unmatched2) = $gene_comparison->get_unmatched_Genes;
#
## this returns an array of GeneCluster objects as well, but only containing the unmatched ones



## get the list of fragmented genes
#
# my @fragmented = $gene_comparison->get_fragmented_Genes (@clusters);
#
## this returns an array of GeneCluster objects as well, but only containing the fragmented ones



## cluster the transcripts using the gene clusters obtained above:
#
# my @transcript_clusters = $gene_comparison->cluster_Transcripts_by_Gene(@clusters);
#
## this returns an array of TranscriptCluster objects, 
## which can be printed as we did with the GeneCluster objects above



###cluster the transcripts of the genes in _gene_array1 and gene_array2 directly
# 
# my @same_transcript_clusters = $gene_comparison->cluster_Transcripts;
#
## this returns an array of TranscriptCluster objects as well






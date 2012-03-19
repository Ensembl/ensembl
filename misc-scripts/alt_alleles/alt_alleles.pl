#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config pass_through);

# (make sure api version is correct
# Usage:
# perl alt_alleles.pl -cpass XXXX > & human_release_63_alt_alleles
#
#
# long way
# perl alt_alleles.pl -vhost ens-staging1 -vport 3306 -vdbname homo_sapiens_vega_63_37 -cdbname homo_sapiens_core_63_37 -chost ens-staging1 -cpass XXXX > & human_release_63_alt_alleles
#

my ($vhost, $vpass, $vport, $vdbname, $vuser, $chost, $cpass, $cport, $cdbname, $cuser);

GetOptions(
    'vuser=s'  => \$vuser,
    'vpass=s'  => \$vpass,
    'vhost=s'  => \$vhost,
    'vport=i'  => \$vport,
    'vdbname=s'       => \$vdbname,
    'cuser=s'  => \$cuser,
    'cpass=s'  => \$cpass,
    'chost=s'  => \$chost,
    'cport=i'  => \$cport,
    'cdbname=s'       => \$cdbname);
#
# Connect to the vgea databse to get the alt allele data.
#

my $api_version = Bio::EnsEMBL::ApiVersion->software_version();

if(!defined($vdbname)){
  $vdbname = "homo_sapiens_vega_".$api_version."_37";
}

if(!defined($cdbname)){
  $cdbname = "homo_sapiens_core_".$api_version."_37";
}

#
# Connect to the core database 
#

my $core_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host => $chost||'ens-staging1',
                                          -user => $cuser||'ensadmin',
					  -pass => $cpass,
                                          -species => "test",
                                          -dbname => $cdbname||"homo_sapiens_core_63_37");



#
# get ensembl gene ids and vega stable ids from the core database
# 


my %vega_to_ens_id;
my ($vega_stable_id, $gene_id);

my $sth =  $core_dba->dbc->prepare("select ensembl_id, display_label from object_xref join xref using(xref_id) join external_db using(external_db_id) where db_name = 'OTTG' and ensembl_object_type = 'Gene'");
$sth->execute;

$sth->bind_columns(\$gene_id, \$vega_stable_id);

while ($sth->fetch){
  $vega_to_ens_id{$vega_stable_id} = $gene_id;
}
$sth->finish;


my $vega_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host => $vhost||'ens-staging1',
                                          -user => $vuser||'ensro',
                                          -port => $vport||3306,
                                          -dbname => $vdbname||"homo_sapiens_vega_63_37");

#
# SQL to get alt_allele data from vega
#

my $sql =(<<EOS);
select aa.alt_allele_id, g.stable_id
 from alt_allele aa, gene g 
 where aa.gene_id = g.gene_id
EOS


#
# Store data in a hash where the key is the alt_id and the ensembl gene ids 
# stored in an anonymous array (value of the hash).
#

my $sth = $vega_dba->dbc->prepare($sql);
$sth->execute;
my ($alt_id, $vega_stable_id);
my %alt_alleles;
$sth->bind_columns(\$alt_id, \$vega_stable_id);


my %no_gene_id;

while($sth->fetch()){
    my $gene_id = $vega_to_ens_id{$vega_stable_id};
  if ( defined($gene_id) ) {
      push @{$alt_alleles{$alt_id}}, $gene_id ;
  } else {
      push @{$no_gene_id{$alt_id}}, $vega_stable_id;
      print STDERR "no ensembl gene_id found for vega stable id $vega_stable_id in core\n";
  }

}
$sth->finish;


#
# Delete the old data
#

my $sth = $core_dba->dbc->prepare("delete from alt_allele");
$sth->execute;


#
# Store alt_alleles.
#

my $alt_allele_count=0;
my $gene_count = 0;

my $ga = $core_dba->get_adaptor("gene");

foreach my $key (keys %alt_alleles){
  my @gene_ids = @{$alt_alleles{$key}};
  my @genes;
  foreach my $gene_id (@gene_ids) {
      push @genes, $ga->fetch_by_dbID($gene_id);
  }

  my $alt_allele_id = $ga->store_alt_alleles(\@genes);
  $alt_allele_count ++ if ($alt_allele_id);
  $gene_count += scalar(@genes) if ($alt_allele_id);
}

print "Added $alt_allele_count alt_allele ids for $gene_count genes\n";

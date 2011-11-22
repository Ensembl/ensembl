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
my %alt_to_gene;
$sth->bind_columns(\$alt_id, \$vega_stable_id);

my $max_alt_id = 0;

my %no_gene_id;

while($sth->fetch()){
  my $gene_id = $vega_to_ens_id{$vega_stable_id};
  if ( defined($gene_id) ) {
      push @{$alt_to_gene{$alt_id}}, $gene_id ;
  } else {
      push @{$no_gene_id{$alt_id}}, $vega_stable_id;
      print STDERR "no ensembl gene_id found for vega stable id $vega_stable_id in core\n";
  }
  if($alt_id > $max_alt_id){
    $max_alt_id = $alt_id;
  }
}
$sth->finish;

$max_alt_id += 1000;


my %non_reference; # $non_reference{$gene_id} = 1;


$sql = (<<SEQ);
SELECT g.gene_id 
  FROM gene g, seq_region_attrib sra, attrib_type at
    WHERE g.seq_region_id = sra.seq_region_id AND
          at.attrib_type_id = sra.attrib_type_id AND
          at.code = 'non_ref'
SEQ

$sth = $core_dba->dbc->prepare($sql);
$sth->execute;
$sth->bind_columns(\$gene_id);
while($sth->fetch()){
  $non_reference{$gene_id} = 1;
}


#
# Delete the old data
#

$sth = $core_dba->dbc->prepare("delete from alt_allele");
$sth->execute;


#
# Check that all alt_alleles are valid.
# NOTE: if an alt allele has only one gene_id delete it from the hash 
#

foreach my $key (keys %alt_to_gene){
  my @arr = @{$alt_to_gene{$key}};
  if(scalar(@arr) <= 1){
      delete $alt_to_gene{$key};
      my @unmapped_vega = $no_gene_id{$key};
      print STDERR "ignoring alt allele $key: unmapped vega stable id(s): ". join(", ",@unmapped_vega) ."\n";    
  }
}

#
# Then store the alt_allele and gene_id in the alt_allele table.
#

$sql = "insert into alt_allele (alt_allele_id, gene_id, is_ref) values(?, ?, ?)";
my $insert_sth = $core_dba->dbc->prepare($sql);

my $count=0;
my $alt_id_count = 0;


#
# gene_id_to_alt_id needed for LRG stuff
#
#my %gene_id_to_alt_id;

foreach my $key (keys %alt_to_gene){
  my @arr = @{$alt_to_gene{$key}};
  if(scalar(@arr) > 1){
    my $sanity_check =0; # should be 1 refernce gene only if 0 or more than 1 we may have a problem
    foreach my $gene_id (@arr){
      my $ref = 1;
      if(defined( $non_reference{$gene_id})){
	$ref = 0;
      }      
      else{
	$sanity_check++;
      }
      $insert_sth->execute($key, $gene_id, $ref);
#      $gene_id_to_alt_id{$gene_id}= $key;
      $count++;
    }
    if($sanity_check !=1){
      print STDERR "Problem we have $sanity_check genes on the reference sequence for alt_id $key\n";
    }
    $alt_id_count++;
  }
  elsif(scalar(@arr) == 0){ # been removed due to stable id error
  }
}


print "Added $alt_id_count alt_allele ids for $count genes\n";

#
# useful sql to look at sanity check problems
#
#select g.stable_id,g.seq_region_start, g.seq_region_end, s.name, x.display_label from gene g, seq_region s, alt_allele aa, xref x where g.display_xref_id = x.xref_id and g.gene_id = aa.gene_id and g.seq_region_id = s.seq_region_id and aa.alt_allele_id =37;



exit;









## LRG SQL. How to fit this in?
##select ox.ensembl_id, g.gene_id from xref x, object_xref ox, external_db e, gene g where x.xref_id = ox.xref_id and e.external_db_id = x.external_db_id and e.db_name like "Ens_Hs_gene" and ox.ensembl_object_type = "Gene" and x.display_label = g.stable_id ;



##
## Use $max_alt_id for new ones.
##

#$sql =(<<LRG);
#SELECT  ox.ensembl_id, g.gene_id 
#  FROM xref x, object_xref ox, external_db e, gene g
#    WHERE x.xref_id = ox.xref_id AND
#          e.external_db_id = x.external_db_id AND
#          e.db_name like "Ens_Hs_gene" AND
#          ox.ensembl_object_type = "Gene" AND
#           x.display_label = g.stable_id
#LRG

#$sth = $core_dba->dbc->prepare($sql);
#my ($core_gene_id, $lrg_gene_id);
#$sth->execute();
#$sth->bind_columns(\$lrg_gene_id, \$core_gene_id);

#$count =0;

#my $old_count = 0;
#my $new_count = 0;
#my $lrg_count = 0;
##
## If the core gene is already in an alt_allele set then use that alt_id for the LRG gene only.
## Else use a new one and add both core and LRG.
##


#while ($sth->fetch()){
#  if(defined($gene_id_to_alt_id{$core_gene_id})){
#    $insert_sth->execute($gene_id_to_alt_id{$core_gene_id}, $lrg_gene_id);
#    $old_count++;
#  }
#  elsif(defined($gene_id_to_alt_id{$lrg_gene_id})){
#    $insert_sth->execute($gene_id_to_alt_id{$lrg_gene_id}, $core_gene_id);
#    print "LRG perculiarity\t$core_gene_id\t$lrg_gene_id\n";
#    $lrg_count++;
#  }
#  else{ # new one.
#    $max_alt_id++;
#    $insert_sth->execute($max_alt_id, $lrg_gene_id);
#    $insert_sth->execute($max_alt_id, $core_gene_id);
#    $new_count++;
#  }
#  $count++;
#}

#print "Added $count alt_allels for the lrgs. $old_count added to previous alt_alleles and $new_count new ones\n";
#print "LRG problem count = $lrg_count\n";

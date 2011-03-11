use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config pass_through);

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
my $vega_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host => $vhost||'ensdb-1-11',
                                          -user => $vuser||'ensro',
                                          -port => $vport||5317,
                                          -dbname => $vdbname||"vega_homo_sapiens_20100903_v61_GRCh37");



#
# the magic SQL to get the data from vega
#

my $sql =(<<EOS);
select aa.alt_allele_id, gsi.stable_id, x2.display_label
 from alt_allele aa, gene_stable_id gsi, xref x1, gene g
     left join object_xref ox on g.gene_id = ox.ensembl_id
     left join xref x2 on ox.xref_id = x2.xref_id
     left join external_db edb on x2.external_db_id = edb.external_db_id
 where aa.gene_id = gsi.gene_id
  and gsi.gene_id = g.gene_id
  and g.display_xref_id = x1.xref_id
  and ox.ensembl_object_type = 'Gene'
  and edb.db_name = 'ENSG';
EOS

#
# Store data in a hash where the key is the alt_id and the stable ids 
# stored in an anonymous array (value of the hash).
#

my $sth = $vega_dba->dbc->prepare($sql);
$sth->execute;
my ($alt_id, $vega_stable_id, $core_stable_id);
my %alt_to_stable;
$sth->bind_columns(\$alt_id, \$vega_stable_id, \$core_stable_id);

while($sth->fetch()){
  push @{$alt_to_stable{$alt_id}}, $core_stable_id;
}
$sth->finish;


#
# Test case for invalid stable id. (Delete after testing
#

push @{$alt_to_stable{9999}}, "INVALID1";

# initial testing to make sure sensible info was obtained.
#foreach my $key (keys %alt_to_stable){
#  my @arr = @{$alt_to_stable{$key}};
#  if(scalar(@arr) > 1){
#    print $key;
#    foreach my $stable_id (@arr){
#      print "\t".$stable_id;
#    }
#    print "\n";
#  }
#  else{
#    print "$key has only one ensembl stable id ".$arr[0]." so ignoring\n";
#  }
#}


#EXAMPLE:
#54      ENSG00000206285 ENSG00000236802 ENSG00000226936 ENSG00000235155 ENSG00000235863


#
# Connect to the core database to store the data in.
#

my $core_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host => $chost||'ens-research',
                                          -user => $cuser||'ensadmin',
					  -pass => $cpass,
                                          -species => "test",
                                          -dbname => $cdbname||"ianl_homo_sapiens_core_61_37f");




#
# for each stable_id look it up in the gene_stable_id table and get the gene_id
# store this in a hash. 

my %stable_id_to_gene_id;

$sth =  $core_dba->dbc->prepare("Select stable_id, gene_id from gene_stable_id");
$sth->execute;
my ($gene_id);
$sth->bind_columns(\$core_stable_id, \$gene_id);

while ($sth->fetch){
  $stable_id_to_gene_id{$core_stable_id} = $gene_id;
}
$sth->finish;


#
# Delete the old data
#

$sth = $core_dba->dbc->prepare("delete from alt_allele");
$sth->execute;


#
# Check that all stable_ids are valid.
# NOTE: if one is invalid then if the alt_allele only had two members
# then with the reduction of one member will invalidate the alt_allele.
#

foreach my $key (keys %alt_to_stable){
  my @arr = @{$alt_to_stable{$key}};
  my @arr2;
  if(scalar(@arr) > 1){
    foreach my $stable_id (@arr){
      if(defined($stable_id_to_gene_id{$stable_id})){
	push @arr2, $stable_id;
      }
      else{
	print STDERR "$stable_id in vega database but not in core\n";
      }
    }
  }
  else{
    print "$key has only one ensembl stable id ".join(", ",@arr)." so ignoring\n";
  }
    $alt_to_stable{$key} = \@arr2;
}

#
# Then store the alt_allele and gene_id in the alt_allele table.
#

$sql = "insert into alt_allele (alt_allele_id, gene_id) values(?, ?)";
$sth = $core_dba->dbc->prepare($sql);

my $count=0;
my $alt_id_count = 0;
foreach my $key (keys %alt_to_stable){
  my @arr = @{$alt_to_stable{$key}};
  if(scalar(@arr) > 1){
    foreach my $stable_id (@arr){
      $sth->execute($key, $stable_id_to_gene_id{$stable_id});
      $count++;
    }
    $alt_id_count++;
  }
  elsif(scalar(@arr) == 0){ # been removed due to stable id error
  }
}


print "Added $alt_id_count alt_allele ids for $count genes\n";





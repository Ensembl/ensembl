#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::AltAlleleGroup;
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
# Connect to the core & vega database 
#

my $core_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host => $chost||'ens-staging1',
  -user => $cuser||'ensadmin',
  -pass => $cpass,
  -group => 'core',
  -dbname => $cdbname,
  -port => $cport
);

my $vega_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host => $vhost||'ens-staging1',
  -user => $vuser||'ensadmin',
  -pass => $vpass,
  -group => 'vega',
  -dbname => $vdbname,
  -port => $vport
);


#
# get ensembl gene ids and vega stable ids from the *core* database
# 
my $vega_core_sql = <<'SQL';
select display_label, ensembl_id
from object_xref 
join xref using(xref_id) 
join external_db using(external_db_id) 
where db_name = 'OTTG' 
and ensembl_object_type = 'Gene'
SQL

# sometimes we will see more than one gene associated with an OTTG
# this happens when an OTTG on the primary assemby has been projected to a patch.
my %vega_to_ensembl_core_gene_id;
$core_dba->dbc->sql_helper()->execute_no_return(-SQL => $vega_core_sql, -CALLBACK => sub {
  my ($row) = @_;
  my ($vega_stable_id, $gene_id) = @{$row};
  $vega_to_ensembl_core_gene_id{$vega_stable_id}{$gene_id} = $gene_id;
});

print "\nFetched ".(scalar(keys %vega_to_ensembl_core_gene_id))." Vega Stable IDs\n";

#
# Get AltAlleles from vega
#
my $vega_aaga = $vega_dba->get_AltAlleleGroupAdaptor();
my $vega_groups = $vega_aaga->fetch_all();

my $cnt_vega_rows = @{$vega_groups};
print STDERR "Fetched $cnt_vega_rows rows from the vega db alt_allele table\n";

my %no_gene_id;
my @new_groups;
foreach my $group (@{$vega_groups}) {
  my $members = $group->get_all_Genes_types();
  my $new_core_group = Bio::EnsEMBL::AltAlleleGroup->new();
  foreach my $member (@{$members}) {
    my ($vega_gene, $attribs_hash) = @{$member};
    my $vega_stable_id = $vega_gene->stable_id();
    if(exists $vega_to_ensembl_core_gene_id{$vega_stable_id}) {
      foreach my $gene_id (keys %{$vega_to_ensembl_core_gene_id{$vega_stable_id}} ) {
        #Add each gene in. If we had a 1:m relationship then we copy the attribute already assigned
        #across
        $new_core_group->add_member($gene_id, $attribs_hash);
      }
    }
    else {
      push @{$no_gene_id{$group->dbID()}}, $vega_stable_id;
      print STDERR "no ensembl gene_id found for vega stable id $vega_stable_id in core\n";
    }
  }
  if($new_core_group->size() > 0) {
    push(@new_groups, $new_core_group);
  }
}

#
# Delete the old data
#
print STDERR "\n\nDeleting all alt_alleles...\n\n";
$core_dba->dbc->do("delete from alt_allele");
$core_dba->dbc->do("delete from alt_allele_attrib");
$core_dba->dbc->do("delete from alt_allele_group");

#
# Store alt_alleles.
#
print STDERR "Storing new alt alleles...\n\n";
my $alt_allele_count=0;
my $gene_count = 0;

my $core_aaga = $core_dba->get_AltAlleleGroupAdaptor();
foreach my $group (@new_groups) {
  my $alt_allele_id = $core_aaga->store($group);
  $alt_allele_count++;
  $gene_count += $group->size()
}

print "Added $alt_allele_count alt_allele ids for $gene_count genes\nDONE\n";

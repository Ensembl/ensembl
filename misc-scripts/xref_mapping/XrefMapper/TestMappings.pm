=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefMapper::TestMappings;
use strict;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use warnings;
use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;
use POSIX;

##########################################
# Testing  (may be moved to healthchecks)
##########################################

##### Unlinked entries ##############

# ERRORS
#    dependent_xref            and xref
#    primary_xref              and xref
#    transcript_direct_xref    and xref
#    translation_direct_xref   and xref
#    gene_direct_xref          and xref
#    synonym                   and xref

#    identity_xref             and object_xref
#    go_xref                   and object_xref

#    gene_transcript_translation   and gene_stable_id
#    gene_transcript_translation   and transcript_stable_id
#    gene_transcript_translation   and translation_stable_id

# WARNNGS
#    gene_direct_xref              and gene_stable_id
#    transcript
#    translation

# All object_xref of type go have a go_xref entry


##### Numbers between xref and core (xref and object_xref) are similar

##### if human or mouse check the number of gene name changes.

sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  if (defined($mapper->previous_core)) {
      $self->core($mapper->previous_core);
  } else {
      $self->core($mapper->core);
  }
  $self->xref($mapper->xref);
  $self->mapper($mapper);
  return $self;
}


sub mapper{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_mapper} = $arg );
  return $self->{_mapper};
}


sub unlinked_entries{
  my ($self) = @_;
  
  my $failed = 0;

  my $xref_id;
  my $count;
  my $dbi = $self->xref->dbc;

  my $sth_stat = $dbi->prepare("insert into process_status (status, date) values('tests_started',now())");
  $sth_stat->execute();

  #    dependent_xref            and xref
  my $count_sql = "select count(1) from dependent_xref d left join xref x on d.master_xref_id = x.xref_id where x.xref_id is null";
  my $sql = "select distinct(d.master_xref_id) from dependent_xref d left join xref x on d.master_xref_id = x.xref_id where x.xref_id is null limit 10";
  my $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print STDERR "SQL QUERY: $sql\n";
    while($sth->fetch){
      print STDERR "Problem with master xref $xref_id\n";
    }
    $sth->finish;
  }




  $count_sql = "select count(1) from dependent_xref d left join xref x on d.dependent_xref_id = x.xref_id where x.xref_id is null";

  $sql = "select distinct(d.dependent_xref_id) from dependent_xref d left join xref x on d.dependent_xref_id = x.xref_id where x.xref_id is null limit 10";

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print STDERR "SQL QUERY: $sql\n";
    while($sth->fetch){
      print STDERR "Problem with dependent xref $xref_id\n";
    }
    $sth->finish;
  }

  $count_sql = "select count(1) from primary_xref d left join xref x on d.xref_id = x.xref_id where x.xref_id is null";

  $sql = "select distinct(d.xref_id) from primary_xref d left join xref x on d.xref_id = x.xref_id where x.xref_id is null limit 10";

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print STDERR "SQL QUERY: $sql\n";
    while($sth->fetch){
      print STDERR "Problem with primary xref $xref_id\n";
    }
    $sth->finish;
  }

  foreach my $type (qw(transcript translation gene)){
    $count_sql =      "select count(1) from ".$type."_direct_xref d left join xref x on d.general_xref_id = x.xref_id where x.xref_id is null";
    
    $sql = "select distinct(d.general_xref_id) from ".$type."_direct_xref d left join xref x on d.general_xref_id = x.xref_id where x.xref_id is null limit 10";
    
    $sth = $dbi->prepare($count_sql);
    $sth->execute();
    $sth->bind_columns(\$count);
    $sth->fetch();
    $sth->finish;
    
    if($count){
      $failed = 1;
      $sth = $dbi->prepare($sql);
      $sth->execute();
      $sth->bind_columns(\$xref_id);
      print STDERR "SQL QUERY: $sql\n";
      while($sth->fetch){
	print STDERR "Problem with ".$type."_direct_xref $xref_id\n";
      }
      $sth->finish;
    }
    
  }


  $count_sql = "select count(1) from synonym d left join xref x on d.xref_id = x.xref_id where x.xref_id is null";

  $sql = "select distinct(d.xref_id) from synonym d left join xref x on d.xref_id = x.xref_id where x.xref_id is null limit 10";

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print STDERR "SQL QUERY: $sql\n";
    while($sth->fetch){
      print STDERR "Problem with synonym $xref_id\n";
    }
    $sth->finish;
  }


  $count_sql = "select count(1) from identity_xref d left join object_xref o on d.object_xref_id = o.object_xref_id where o.object_xref_id is null";

  $sql = "select distinct(d.object_xref_id) from identity_xref d left join object_xref o on d.object_xref_id = o.object_xref_id where o.object_xref_id is null limit 10";

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print STDERR "SQL QUERY: $sql\n";
    while($sth->fetch){
      print STDERR "Problem with object_xref $xref_id\n";
    }
    $sth->finish;
  }

  $count_sql = "select count(1) from go_xref d left join object_xref o on d.object_xref_id = o.object_xref_id where o.object_xref_id is null";

  $sql = "select distinct(d.object_xref_id) from go_xref d left join object_xref o on d.object_xref_id = o.object_xref_id where o.object_xref_id is null limit 10";

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print STDERR "SQL QUERY: $sql\n";
    while($sth->fetch){
      print STDERR "Problem with object_xref $xref_id\n";
    }
    $sth->finish;
  }


  foreach my $type (qw(transcript translation gene)){
    $count_sql = "select count(1) from gene_transcript_translation d left join ".$type."_stable_id x on d.".$type."_id = x.internal_id where x.internal_id is null and d.".$type."_id is not null";
    
    $sql = "select distinct(d.".$type."_id) from gene_transcript_translation d left join  ".$type."_stable_id x on d.".$type."_id = x.internal_id where x.internal_id is null and d.".$type."_id is not null limit 10";
    
    $sth = $dbi->prepare($count_sql);
    $sth->execute();
    $sth->bind_columns(\$count);
    $sth->fetch();
    $sth->finish;
    
    if($count){
      $failed = 1;
      $sth = $dbi->prepare($sql);
      $sth->execute();
      $sth->bind_columns(\$xref_id);
      print STDERR "SQL QUERY: $sql\n";
      while($sth->fetch){
	print STDERR "Problem with ".$type."_id $xref_id\n";
      }
      $sth->finish;
    }
    
  }


  $count_sql = "select count(1) from xref x, source s, object_xref o left join go_xref g on o.object_xref_id = g.object_xref_id where x.xref_id = o.xref_id and s.source_id = x.source_id and s.name like 'GO' and ox_status in ('DUMP_OUT') and g.object_xref_id is null";
  $sql = "select distinct(o.object_xref_id) from xref x, source s, object_xref o left join go_xref g on o.object_xref_id = g.object_xref_id where x.xref_id = o.xref_id and s.source_id = x.source_id and s.name like 'GO' and ox_status in ('DUMP_OUT') and g.object_xref_id is null limit 10";

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print STDERR "SQL QUERY: $sql\n";
    while($sth->fetch){
      print STDERR "Problem with object_xref $xref_id which is linked to a GO source but has no go_xref reference\n";
    }
    $sth->finish;
  }

  if(!$failed){
    $sth_stat = $dbi->prepare("insert into process_status (status, date) values('tests_finished',now())");
    $sth_stat->execute();
  }
  else{
    $sth_stat = $dbi->prepare("insert into process_status (status, date) values('tests_failed',now())");
    $sth_stat->execute();
  }
  $sth_stat->finish;

  return $failed;
}


sub entry_number_check{
  my ($self) = @_;
  

# No point doing xrefs object_xrefs are more important and gives a better indication of wether things went okay.

  my %old_object_xref_count;
  my %new_object_xref_count;
  my $dbi = $self->xref->dbc;

  my $sth = $dbi->prepare('select s.name, count(distinct x.xref_id, ensembl_id) from xref x, object_xref ox, source s where ox.xref_id = x.xref_id  and x.source_id = s.source_id and ox_status = "DUMP_OUT" and s.name not like "AFFY%"   group by s.name');
  $sth->execute();
  my ($name, $count);
  $sth->bind_columns(\$name,\$count);
  while($sth->fetch()){
    $new_object_xref_count{$name} = $count;
  }
  $sth->finish;

  
 $sth = $self->core->dbc->prepare('select e.db_name, count(*) from xref x, object_xref ox, external_db e where ox.xref_id = x.xref_id and x.external_db_id = e.external_db_id and e.db_name not like "AFFY%" and (x.info_type is NULL or x.info_type != "PROJECTION") group by e.db_name');

  $sth->execute();
  $sth->bind_columns(\$name,\$count);
  while($sth->fetch()){
    my $change = 0;
    $old_object_xref_count{$name} = $count;
    if(defined($new_object_xref_count{$name})){
      $change = (($new_object_xref_count{$name} - $count)/$count) * 100;
      if($change > 5){ # increase of 5%
	print "WARNING: $name has increased by ".int($change)."\% was $count now ". $new_object_xref_count{$name} . "\n" if($self->mapper->verbose); 
      }
      elsif($change < -5){ # decrease by 5%
	print "WARNING: $name has decreased by ".int($change)." \% was $count now ". $new_object_xref_count{$name} . "\n" if($self->mapper->verbose); 
      }
    }
    else{
      print "WARNING: xrefs $name are not in the new database but $count are in the old???\n" if($self->mapper->verbose);
    }
  }
  $sth->finish;
  
  foreach my $key (keys %new_object_xref_count){
    if(!defined($old_object_xref_count{$key})){
      print "WARNING: $key has ".$new_object_xref_count{$key} ." xrefs in the new database but NONE in the old\n" if($self->mapper->verbose);
    }
  }
  
  return;
}


sub name_change_check{
  my ($self) = @_;

  my %new_name; # $old_name{$gene_id} = HGNC_%name
  my %id_to_stable_id;
  my $dbi = $self->xref->dbc;

  my $official_name = $self->mapper->get_official_name;
  if(!defined($official_name)){
    return;
  }
#  print "Checking names\n";

  my $sql = 'select x.label, gsi.internal_id, gsi.stable_id from object_xref ox, xref x, gene_stable_id gsi, source s  where x.xref_id = ox.xref_id and ox.ensembl_object_type = "Gene" and gsi.internal_id = ox.ensembl_id and x.source_id = s.source_id and s.name like "'.$official_name.'_%"';

  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($name, $gene_id, $stable_id);
  $sth->bind_columns(\$name,\$gene_id, \$stable_id);
  my $count = 0;
  while($sth->fetch()){
    $new_name{$gene_id} = $name;
    $id_to_stable_id{$gene_id} = $stable_id;
    $count++;
  }
  $sth->finish;
#  print $sql."\n";
#  print $count." entries found in xref database\n";



  # Use synonyms as well.
  my %alias;
  $sql = 'select x.label, sy.synonym from xref x, synonym sy, source so where x.xref_id = sy.xref_id and so.source_id = x.source_id and so.name like "'.$official_name.'_%" ';
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($syn);
  $sth->bind_columns(\$name,\$syn);
  $count = 0;
  while($sth->fetch()){
    $alias{$syn} = $name;
  }
  $sth->finish;  

  $sql = 'select x.label, sy.synonym from xref x, synonym sy, source so where x.xref_id = sy.xref_id and so.source_id = x.source_id and so.name like "EntrezGene"';
  $sth = $dbi->prepare($sql);
  $sth->execute();
  $sth->bind_columns(\$name,\$syn);
  while($sth->fetch()){
    $alias{$syn} = $name;
  }
  $sth->finish;



  # NOTE ncRNA has higher priority
  $sql = "select x.display_label, g.gene_id from gene g, xref x where g.display_xref_id = x.xref_id and biotype = 'protein_coding'";

  $sth = $self->core->dbc->prepare($sql);
  $sth->execute();
  $sth->bind_columns(\$name,\$gene_id);
  $count =0;
  my $total_count=0;
  while($sth->fetch()){
    if(defined($new_name{$gene_id})){
      $total_count++;
    }	
    if(defined($new_name{$gene_id}) and $new_name{$gene_id} ne $name){
      if(!defined($alias{$name}) or $alias{$name} ne $new_name{$gene_id}){
	print STDERR "WARN: gene_id ($gene_id) ".$id_to_stable_id{$gene_id}." new = ".$new_name{$gene_id}." old = $name\n";
	$count++;
      }
    }	
  }
  if($total_count){
    print STDERR "$count entries with different names out of $total_count protein coding gene comparisons?\n";
  }
}


sub direct_stable_id_check{
  my ($self) = @_;

  my $dbi = $self->xref->dbc;

  foreach my $type (qw(gene transcript translation)){
    
    my $sql = "select s.name, count(*) from source s, xref x, ".$type."_direct_xref gdx left join ".$type."_stable_id gsi on gdx.ensembl_stable_id = gsi.stable_id where s.source_id = x.source_id and x.xref_id = gdx.general_xref_id and gsi.stable_id is null group by s.name";
    
    my $sth = $dbi->prepare($sql);
    $sth->execute();
    my ($name, $count);
    $sth->bind_columns(\$name,\$count);
    my $total_count=0;
    while($sth->fetch()){
      print STDERR "WARNING $name has $count invalid stable ids in ".$type."_direct_xrefs\n";
      $total_count += $count;
    }	
    $sth->finish;
    if($total_count){
      print STDERR "USEFUL SQL: $sql\n";
    }
  }

}
1;

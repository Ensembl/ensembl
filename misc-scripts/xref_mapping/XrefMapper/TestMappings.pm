package XrefMapper::TestMappings;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use strict;
use warnings;
use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;

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

# All object_xref of type go have a go_xref entry


##### Numbers between xref and core (xref and object_xref) are similar

##### if human or mouse check the number of gene name changes.

sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
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

  #    dependent_xref            and xref
  my $count_sql = "select count(1) from dependent_xref d left join xref x on d.master_xref_id = x.xref_id where x.xref_id is null";
  my $sql = "select distinct(d.master_xref_id) from dependent_xref d left join xref x on d.master_xref_id = x.xref_id where x.xref_id is null limit 10";
  my $sth = $self->xref->dbc->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print "SQL QUERY: $sql\n";
    while($sth->fetch){
      print "Problem with master xref $xref_id\n";
    }
    $sth->finish;
  }




  $count_sql = "select count(1) from dependent_xref d left join xref x on d.dependent_xref_id = x.xref_id where x.xref_id is null";

  $sql = "select distinct(d.dependent_xref_id) from dependent_xref d left join xref x on d.dependent_xref_id = x.xref_id where x.xref_id is null limit 10";

  $sth = $self->xref->dbc->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print "SQL QUERY: $sql\n";
    while($sth->fetch){
      print "Problem with dependent xref $xref_id\n";
    }
    $sth->finish;
  }

  $count_sql = "select count(1) from primary_xref d left join xref x on d.xref_id = x.xref_id where x.xref_id is null";

  $sql = "select distinct(d.xref_id) from primary_xref d left join xref x on d.xref_id = x.xref_id where x.xref_id is null limit 10";

  $sth = $self->xref->dbc->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print "SQL QUERY: $sql\n";
    while($sth->fetch){
      print "Problem with primary xref $xref_id\n";
    }
    $sth->finish;
  }

  foreach my $type (qw(transcript translation gene)){
    $count_sql = "select count(1) from ".$type."_direct_xref d left join xref x on d.general_xref_id = x.xref_id where x.xref_id is null";
    
    $sql = "select distinct(d.xref_id) from ".$type."_direct_xref d left join xref x on d.general_xref_id = x.xref_id where x.xref_id is null limit 10";
    
    $sth = $self->xref->dbc->prepare($count_sql);
    $sth->execute();
    $sth->bind_columns(\$count);
    $sth->fetch();
    $sth->finish;
    
    if($count){
      $failed = 1;
      $sth = $self->xref->dbc->prepare($sql);
      $sth->execute();
      $sth->bind_columns(\$xref_id);
      print "SQL QUERY: $sql\n";
      while($sth->fetch){
	print "Problem with ".$type."_direct_xref $xref_id\n";
      }
      $sth->finish;
    }
    
  }


  $count_sql = "select count(1) from synonym d left join xref x on d.xref_id = x.xref_id where x.xref_id is null";

  $sql = "select distinct(d.xref_id) from synonym d left join xref x on d.xref_id = x.xref_id where x.xref_id is null limit 10";

  $sth = $self->xref->dbc->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print "SQL QUERY: $sql\n";
    while($sth->fetch){
      print "Problem with synonym $xref_id\n";
    }
    $sth->finish;
  }


  $count_sql = "select count(1) from identity_xref d left join object_xref o on d.object_xref_id = o.object_xref_id where o.object_xref_id is null";

  $sql = "select distinct(d.object_xref_id) from identity d left join object_xref o on d.object_xref_id = o.object_xref_id where o.object_xref_id is null limit 10";

  $sth = $self->xref->dbc->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print "SQL QUERY: $sql\n";
    while($sth->fetch){
      print "Problem with object_xref $xref_id\n";
    }
    $sth->finish;
  }

  $count_sql = "select count(1) from go_xref d left join object_xref o on d.object_xref_id = o.object_xref_id where o.object_xref_id is null";

  $sql = "select distinct(d.object_xref_id) from go_xref d left join object_xref o on d.object_xref_id = o.object_xref_id where o.object_xref_id is null limit 10";

  $sth = $self->xref->dbc->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print "SQL QUERY: $sql\n";
    while($sth->fetch){
      print "Problem with object_xref $xref_id\n";
    }
    $sth->finish;
  }


  foreach my $type (qw(transcript translation gene)){
    $count_sql = "select count(1) from gene_transcript_translation d left join ".$type."_stable_id x on d.".$type."_id = x.internal_id where x.internal_id is null and d.".$type."_id is not null";
    
    $sql = "select distinct(d.".$type."_id) from gene_transcript_translation d left join  ".$type."_stable_id x on d.".$type."_id = x.internal_id where x.internal_id is null and d.".$type."_id is not null limit 10";
    
    $sth = $self->xref->dbc->prepare($count_sql);
    $sth->execute();
    $sth->bind_columns(\$count);
    $sth->fetch();
    $sth->finish;
    
    if($count){
      $failed = 1;
      $sth = $self->xref->dbc->prepare($sql);
      $sth->execute();
      $sth->bind_columns(\$xref_id);
      print "SQL QUERY: $sql\n";
      while($sth->fetch){
	print "Problem with ".$type."_id $xref_id\n";
      }
      $sth->finish;
    }
    
  }


  $count_sql = "select count(1) from xref x, source s, object_xref o left join go_xref g on o.object_xref_id = g.object_xref_id where x.xref_id = o.xref_id and s.source_id = x.source_id and s.name like 'GO' and g.object_xref_id is null";
  $sql = "select distinct(o.object_xref_id) from xref x, source s, object_xref o left join go_xref g on o.object_xref_id = g.object_xref_id where x.xref_id = o.xref_id and s.source_id = x.source_id and s.name like 'GO' and g.object_xref_id is null limit 10";
  my $sth = $self->xref->dbc->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns(\$count);
  $sth->fetch();
  $sth->finish;
  
  if($count){
    $failed = 1;
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print "SQL QUERY: $sql\n";
    while($sth->fetch){
      print "Problem with object_xref $xref_id which is linked to a GO source but has no go_xref reference\n";
    }
    $sth->finish;
  }




  return $failed;
}


sub entry_number_check{
  my ($self) = @_;
  

# No point doing xrefs object_xrefs are more important and gives a better indication of wether things went okay.

  my %old_object_xref_count;
  my %new_object_xref_count;

  my $sth = $self->xref->dbc->prepare('select s.name, count(*) from xref x, object_xref ox, source s where ox.xref_id = x.xref_id  and x.source_id = s.source_id and ox_status = "DUMP_OUT" group by s.name');
  $sth->execute();
  my ($name, $count);
  $sth->bind_columns(\$name,\$count);
  while($sth->fetch()){
    $new_object_xref_count{$name} = $count;
  }
  $sth->finish;

  
 $sth = $self->core->dbc->prepare('select e.db_name, count(*) from xref x, object_xref ox, external_db e where ox.xref_id = x.xref_id and x.external_db_id = e.external_db_id group by e.db_name');

  $sth->execute();
  $sth->bind_columns(\$name,\$count);
  while($sth->fetch()){
    my $change = 0;
    if(defined($new_object_xref_count{$name})){
      $change = (($new_object_xref_count{$name} - $count)/$new_object_xref_count{$name}) * 100;
      if($change > 5){ # increase of 5%
	print "WARNING: $name has increased by $change\% was $count now ". $new_object_xref_count{$name}. "\n"; 
      }
      elsif($change < -5){ # decrease by 5%
	print "WARNING: $name has decreased by $change\% was $count now ". $new_object_xref_count{$name}. "\n"; 
      }
    }   
  }
  $sth->finish;

  return;
}


sub name_change_check{
  my ($self) = @_;

  my %old_name; # $old_name{$gene_id} = HGNC_%name


  my $official_name = $self->mapper->get_official_name;
  if(!defined($official_name)){
    return;
  }
  print "Checking names\n";

  my $sql = 'select x.label, gsi.internal_id from object_xref ox, xref x, gene_stable_id gsi, source s  where x.xref_id = ox.xref_id and ox.ensembl_object_type = "Gene" and gsi.internal_id = ox.ensembl_id and x.source_id = s.source_id and s.name like "'.$official_name.'_%"';

  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();
  my ($name, $gene_id);
  $sth->bind_columns(\$name,\$gene_id);
  my $count = 0;
  while($sth->fetch()){
    $old_name{$gene_id} = $name;
  }
  $sth->finish;
  print $sql."\n";
  print $count." etires found in xref database\n";

  $sql = "select x.display_label, g.gene_id from gene g, xref x where g.display_xref_id = x.xref_id";

  $sth = $self->core->dbc->prepare($sql);
  $sth->execute();
  $sth->bind_columns(\$name,\$gene_id);
  $count =0;
  my $total_count=0;
  while($sth->fetch()){
    if(defined($old_name{$gene_id})){
      $total_count++;
    }	
    if(defined($old_name{$gene_id}) and $old_name{$gene_id}ne $name){
      print "WARN: gene_id $gene_id old = ".$old_name{$gene_id}." new = $name\n";
      $count++;
    }	
  }
  print "$count entris with different names out of $total_count?\n";
}


1;

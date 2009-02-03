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

##### if human or mouse check the number of gene name changes.

##### All go_xref have a linkage_type 

##### All object_xref of type go have a go_xref entry

##### Numbers between xref and core (xref and object_xref) are similar


sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  return $self;
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
    $failed = 0;
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    $sth->bind_columns(\$xref_id);
    print "SQL QUERY: $sql\n";
    while($sth->fetch){
      print "Problem with dependent xref $xref_id\n";
    }
    $sth->finish;
  }

  return $failed;
}

1;

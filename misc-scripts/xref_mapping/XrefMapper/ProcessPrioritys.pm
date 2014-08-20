=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package XrefMapper::ProcessPrioritys;
use strict;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use warnings;
use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;

# Process the priority xrefs.

#
# 1) create a list of source "names" that are priority xrefs
#
# 2) Just to be sure set all ox_status in object_xref to 'DUMP_OUT'
#    set dumped in xref to NULL
# 
# 3) for each of the source names 
#    set ox_status to 'FAILED_PRIORITY' for those not the best match
#        Also do this fro its depenedents
#

sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
#  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->verbose($mapper->verbose);
  return $self;
}

sub get_priority_names{
  my ($self) = @_;


  my $psth = $self->xref->dbc->prepare("select s.priority_description, s.name from source s, xref x where x.source_id = s.source_id group by s.priority_description, s.name order by s.name") || die "prepare failed";
  $psth->execute() || die "execute failed";

  my @names;
  my %seen;

  my $last_name = "rubbish";
  my ($desc,$name);
  $psth->bind_columns(\$desc,\$name);
  while($psth->fetch()){
    if($name eq $last_name and !defined($seen{$name})){
      push @names, $name;
      $seen{$name} = 1;
    }
    $last_name = $name;
  }

  return @names;
}


sub process {
  my ($self) = @_;

  my @names = $self->get_priority_names();

  print "The foillowing will be processed as priority xrefs\n" if($self->verbose);
  foreach my $name (@names){
    print "\t$name\n" if($self->verbose);
  }

  my $update_ox_sth = $self->xref->dbc->prepare('update object_xref set ox_status = "FAILED_PRIORITY" where object_xref_id = ?');
  my $update_x_sth  = $self->xref->dbc->prepare("update xref set dumped = 'NO_DUMP_ANOTHER_PRIORITY' where xref_id = ?");
  my $dep_sth       = $self->xref->dbc->prepare("select dependent_xref_id, object_xref_id from dependent_xref where master_xref_id = ?");

  #
  # Change of tact here to make the sql easier...
  #

  # 1) Set to failed all those that have no object xrefs.

  my $f_sql =(<<FSQL);
   SELECT x.xref_id
     FROM  source s, xref x 
     LEFT JOIN object_xref ox ON  ox.xref_id = x.xref_id 
       WHERE x.source_id = s.source_id
         AND s.name = ? 
         AND ox.object_xref_id is null
FSQL

  my $f_sth =  $self->xref->dbc->prepare($f_sql);
  foreach my $name (@names){
    $f_sth->execute($name);
    my ($xref_id);
    $f_sth->bind_columns(\$xref_id);
    while($f_sth->fetch()){
      $update_x_sth->execute($xref_id);
    }
  }
  $f_sth->finish;


  # 
  # Now ALL object_xrefs have an identity_xref :-)
  # So we can do a straight join and treat all info_types the same way.
  # 
  my $new_sql =(<<NEWS);
   SELECT ox.object_xref_id, x.accession, x.xref_id, (ix.query_identity + ix.target_identity) as identity, ox.ox_status
      FROM object_xref ox, xref x, source s, identity_xref ix
        WHERE ox.object_xref_id = ix.object_xref_id 
          AND ox.xref_id = x.xref_id
          AND s.source_id = x.source_id
          AND s.name = ?
         ORDER BY x.accession DESC, s.priority ASC , identity DESC, x.xref_id DESC
NEWS

  my $sth = $self->xref->dbc->prepare($new_sql);
  foreach my $name (@names){
    $sth->execute($name);
    my ($object_xref_id, $acc, $xref_id, $identity, $status);
    $sth->bind_columns(\$object_xref_id, \$acc, \$xref_id, \$identity, \$status);
    my $last_acc = "";
    my $last_name = "";
    my $best_xref_id = undef;
    my @gone;
    while($sth->fetch){
      if($last_acc eq $acc){
	if($xref_id != $best_xref_id){
	  if($status eq "DUMP_OUT"){
	    $update_x_sth->execute($xref_id);
	    $update_ox_sth->execute($object_xref_id);
#??? what about dependents
	    $self->process_dependents($xref_id, $update_ox_sth, $dep_sth);
	  }
	  else{
	    $update_x_sth->execute($xref_id);
	  }
	}
	if(@gone){ #best priority failed so anothre one now found so set dumped;
	  if($last_name eq $acc){
	    foreach my $d (@gone){
	      $update_x_sth->execute($d);
	    }
	  }
	}
      }
      else{ # NEW
	if($status eq "DUMP_OUT"){
	  $last_acc = $acc;
	  $best_xref_id = $xref_id;
	  if(@gone and $last_name eq $acc){
	    foreach my $d (@gone){
	      $update_x_sth->execute($d);
	    }
	    @gone=();
	  }
	}
	else{
	  if($last_name eq $acc){
	    push @gone, $xref_id;	  }
	  else{
	    @gone = ();
	    push @gone, $xref_id;
	  }	
	  $last_name = $acc;
	}
      }
    }
  }
  $sth->finish;

  $update_ox_sth->finish;
  $update_x_sth->finish;


# We want to make sure that if a priority xref is NOT MAPPEd then we only



  $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('prioritys_flagged',now())");
  $sth->execute();
  $sth->finish;
}

sub process_dependents{
  my ($self, $master_xref_id, $update_ox_sth, $dep_sth)= @_;

  my @master_xrefs;

  push @master_xrefs, $master_xref_id;

  while(my $xref_id = pop(@master_xrefs)){
    $dep_sth->execute($xref_id);
    my ($dep_xref_id, $object_xref_id);
    $dep_sth->bind_columns(\$dep_xref_id,\$object_xref_id);
    while($dep_sth->fetch()){
      $update_ox_sth->execute($object_xref_id);
      push @master_xrefs, $dep_xref_id; # remember dependents dependents
    }
  }

}

1;

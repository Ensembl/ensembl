package XrefMapper::ProcessPrioritys;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use strict;
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

  print "The foillowing will be processed as priority xrefs\n";
  foreach my $name (@names){
    print "\t$name\n";
  }

  my $sth =  $self->xref->dbc->prepare('select ox.object_xref_id, x.accession, s.priority, ox.linkage_type from object_xref ox, xref x, source s where x.xref_id = ox.xref_id and x.source_id = s.source_id and s.name = ? and ox.ox_status = "DUMP_OUT" order by x.accession, s.priority');


  my $seq_sth = $self->xref->dbc->prepare('select ox.object_xref_id, (ix.query_identity + ix.target_identity) as identity from object_xref ox, xref x, source s, identity_xref ix where ix.object_xref_id = ox.object_xref_id and x.xref_id = ox.xref_id and ox.ox_status = "DUMP_OUT" and x.source_id = s.source_id and s.name = ? and x.accession = ?  order by identity DESC');

  my $dep_sth = $self->xref->dbc->prepare('select distinct dox.object_xref_id, (ix.query_identity+ ix.target_identity) as identity from dependent_xref dx, object_xref dox, xref x, source s, object_xref mox, identity_xref ix where ix.object_xref_id = mox.object_xref_id and mox.xref_id = dx.master_xref_id and dx.dependent_xref_id = dox.xref_id and dox.xref_id = x.xref_id and s.source_id = x.source_id and x.accession = ? and s.name = ? and dox.ox_status = "DUMP_OUT"');


  my $update_ox_sth = $self->xref->dbc->prepare('update object_xref set ox_status = "FAILED_PRIORITY" where object_xref_id = ?');

  foreach my $name (@names){
    my %priority_clash_seq;
    my %priority_clash_depend;
    print "processing $name source\n";
    $sth->execute($name);
    my ($ox,$acc,$priority, $type);
    $sth->bind_columns(\$ox,\$acc,\$priority,\$type);
    my $last_acc = "";
    my $top_priority = 0;
    while($sth->fetch()){
      if($acc eq $last_acc){
	if($priority == $top_priority and $type ne "DIRECT"){
	  if($type eq "SEQUENCE_MATCH"){
	    $priority_clash_seq{$acc} = 1;
	  }	
	  elsif($type eq "DEPENDENT"){
	    $priority_clash_depend{$acc} = 1;
	  }
	  else{
	    print "type $type????? $acc, $name\n";
	  }
	}
	else{
	  $update_ox_sth->execute($ox);
	}	
      }
      else{
	$top_priority = $priority;
      }
      $last_acc = $acc;
    }	

    # Easy case of sequence matches
    my $identity;
#    print "need to sort out those with the same priority status:-\n";
    foreach my $key (keys %priority_clash_seq){
      $seq_sth->execute($name, $key);
      $seq_sth->bind_columns(\$ox, \$identity);
      $seq_sth->fetch();  # keep first one
      while( $seq_sth->fetch()){
	$update_ox_sth->execute($ox);
      }
    }
    
    # now the hard one dependent xrefs!!
    foreach my $key (keys %priority_clash_depend){
      $dep_sth->execute($key, $name);
      $dep_sth->bind_columns(\$ox, \$identity);
      $dep_sth->fetch();  # keep first one
      while( $dep_sth->fetch()){
	$update_ox_sth->execute($ox); # remove others
      }
    }



  }

  #Sanity check make sure only one instance of each priority xref in object_xref table with ox_status = DUMP_OUT
  foreach my $name (@names){
    print "checking $name source\n";
    $sth->execute($name);
    my ($ox,$acc,$priority, $type);
    $sth->bind_columns(\$ox,\$acc,\$priority,\$type);
    my $last_acc = "";
    while($sth->fetch()){
      if($acc eq $last_acc){
	print "ERROR: $acc still has more than one viable entry in object_xref\n";
      }
      $last_acc = $acc;
    }
  }


  $seq_sth->finish;
  $sth->finish;
  $update_ox_sth->finish;

  $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('prioritys_flagged',now())");
  $sth->execute();
  $sth->finish;
}

1;

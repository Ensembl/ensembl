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
  my $update_x_sth = $self->xref->dbc->prepare('update xref set dumped = 0 where xref_id = ?');

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
   SELECT ox.object_xref_id, x.accession, x.xref_id, (ix.query_identity + ix.target_identity) as identity
      FROM object_xref ox, xref x, source s, identity_xref ix
        WHERE ox.object_xref_id = ix.object_xref_id 
          AND ox.xref_id = x.xref_id
          AND s.source_id = x.source_id
          AND ox.ox_status = "DUMP_OUT"
          AND s.name = ?
         ORDER BY x.accession DESC, s.priority ASC , identity DESC, x.xref_id DESC
NEWS

  my $sth = $self->xref->dbc->prepare($new_sql);
  foreach my $name (@names){
    $sth->execute($name);
    my ($object_xref_id, $acc, $xref_id, $identity);
    $sth->bind_columns(\$object_xref_id, \$acc, \$xref_id, \$identity);
    my $last_acc = "";
    my $best_xref_id = undef;
    while($sth->fetch){
      if($last_acc eq $acc){
         if($xref_id != $best_xref_id){
            $update_x_sth->execute($xref_id);
            $update_ox_sth->execute($object_xref_id);            
         }
      }
      else{ # NEW
        $last_acc = $acc;
        $best_xref_id = $xref_id;
      } 
    }
  }
  $sth->finish;

  $update_ox_sth->finish;
  $update_x_sth->finish;

  $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('prioritys_flagged',now())");
  $sth->execute();
  $sth->finish;
}

1;

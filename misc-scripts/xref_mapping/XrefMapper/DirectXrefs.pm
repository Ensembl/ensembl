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

package XrefMapper::DirectXrefs;
use strict;
use warnings;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;


sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->verbose($mapper->verbose);
  return $self;
}

sub get_ins_ix_sth {
  my $self = shift;
  my $dbi = shift;

  my $sql = (<<"IIX");
INSERT IGNORE INTO identity_xref (object_xref_id, query_identity, target_identity)
  VALUES (?, 100, 100)
IIX
  my $sth = $dbi->prepare($sql);
  return $sth;
}


sub process {
  my $self = shift;

  # Now process the direct xrefs and add data to the object xrefs remember dependent xrefs.

  my $object_xref_id;
  my $dbi = $self->xref->dbc;

  # First get the sths needed for the processing of the direct xrefs;
  my $ins_ox_sql = (<<"IOS");
INSERT INTO object_xref (ensembl_id, xref_id, ensembl_object_type, linkage_type) 
  VALUES (?, ?, ?, ?)
IOS
  my $ins_ox_sth = $dbi->prepare($ins_ox_sql);
  my $get_object_xref_id_sth = $self->get_ox_id_sth($dbi);

  # Direct xrefs can be considered to be 100% matching

  my $ins_ix_sth = $self->get_ins_ix_sth($dbi);

my $stable_sql=(<<"SQL");
  SELECT so.name, dx.general_xref_id, s.internal_id, dx.ensembl_stable_id , dx.linkage_xref
    FROM source so, xref x, TYPE_direct_xref dx left join TYPE_stable_id s on s.stable_id = dx.ensembl_stable_id
      WHERE x.xref_id = dx.general_xref_id and x.source_id = so.source_id 
SQL



  # We want to process the errors ourselves as for greater control
  # If we get a error adding an object xref then it is already there
  # This is not a problem. But we want to know how amny of these there were.

  local $ins_ox_sth->{RaiseError} = 0;  # want to see duplicates and not add de
  local $ins_ox_sth->{PrintError} = 0;

  my %err_count;
 
  foreach my $table (qw(gene transcript translation)){
   my ($dbname, $xref_id, $internal_id, $stable_id, $linkage_type);
   my $sql = $stable_sql;
   $sql =~ s/TYPE/$table/g;
   my $sth = $dbi->prepare($sql);
   $sth->execute();
   $sth->bind_columns(\$dbname, \$xref_id, \$internal_id, \$stable_id, \$linkage_type);
   my $count =0;
   my $duplicate_direct_count = 0;
   my $duplicate_dependent_count = 0;
   while($sth->fetch){
     if(!defined($internal_id)){ # not found either it is an internal id already or stable_id no longer exists
       if($stable_id =~ /^\d+$/){
          $internal_id = $stable_id;
       }
       else{
	 if((!defined($err_count{$dbname})) or ($err_count{$dbname} < 10)){
	   print "Could not find stable id $stable_id in table to get the internal id hence ignoring!!! (for $dbname)\n" if($self->verbose);
	 }
	 $err_count{$dbname}++;
         next;
       }
     }
     $object_xref_id++;
     $count++;
     my @master_xref_ids;
     if($internal_id == 0){
       die "Problem could not find stable id $stable_id and got past the first check for $dbname\n";
     }
     $ins_ox_sth->execute($internal_id, $xref_id, $table, 'DIRECT');
     $get_object_xref_id_sth->execute($internal_id, $xref_id, $table, 'DIRECT');
     $object_xref_id = ($get_object_xref_id_sth->fetchrow_array())[0];
     if($ins_ox_sth->err){
       $duplicate_direct_count++;
       next; #duplicate
     }
     else{
       $ins_ix_sth->execute($object_xref_id);
       push  @master_xref_ids, $xref_id;
     }
     $self->process_dependents({master_xrefs        => \@master_xref_ids,
				dup_count          => \$duplicate_dependent_count,
				table              => $table,
				internal_id        => $internal_id,
                                dbi                => $dbi,
			       });
 
   }
   $sth->finish;
   if($duplicate_direct_count or $duplicate_dependent_count){
     print "duplicate entrys ignored for $duplicate_direct_count direct xrefs and  $duplicate_dependent_count dependent xrefs\n" if($self->verbose);
   }
 }
 foreach my $key ( keys %err_count){
   print STDERR "*WARNING*: ".$err_count{$key}." direct xrefs for database ".$key." could not be added as their stable_ids could not be found\n";
 }

 my $sth = $dbi->prepare("insert into process_status (status, date) values('direct_xrefs_parsed',now())");
 $sth->execute();
 $sth->finish;

 return;
}


sub get_dep_sth {
  my $self = shift;
  my $dbi = shift;

  my $dep_sql = (<<"DSS");
SELECT dependent_xref_id, linkage_annotation
  FROM dependent_xref
    WHERE master_xref_id = ?
DSS
  my $sth = $dbi->prepare($dep_sql);
  return $sth;
}


sub get_dep_go_sth {
  my $self = shift;
  my $dbi = shift;

  my $sql = (<<"IGO");
INSERT INTO go_xref (object_xref_id, linkage_type, source_xref_id) 
  VALUES (?,?,?)
IGO
  my $sth = $dbi->prepare($sql);
  return $sth;
}


sub get_add_dep_ox {
  my $self = shift;
  my $dbi = shift;

  my $sql = (<<"IO2");
INSERT INTO object_xref (ensembl_id, xref_id, ensembl_object_type, linkage_type, master_xref_id)
   VALUES (?, ?, ?, ?, ?)
IO2
  my $sth = $dbi->prepare($sql);
  return $sth;
}

sub get_ox_id_sth {
  my $self = shift;
  my $dbi = shift;

  my $sql = (<<"IO2");
select object_xref_id from object_xref where ensembl_id = ? and xref_id = ? and ensembl_object_type = ? and linkage_type = ?
IO2
  my $sth = $dbi->prepare($sql);
  return $sth;
}

sub get_ox_id_master_sth {
  my $self = shift;
  my $dbi = shift;

  my $sql = (<<"IO2");
select object_xref_id from object_xref where ensembl_id = ? and xref_id = ? and ensembl_object_type = ? and linkage_type = ? and master_xref_id = ?
IO2
  my $sth = $dbi->prepare($sql);
  return $sth;
}


sub process_dependents {
  my ($self, $arg_ref) = @_;

  my $dbi                 = $arg_ref->{dbi};
  my $master_xref_ids     = $arg_ref->{master_xrefs};
  my $duplicate_dep_count = $arg_ref->{dup_count};
  my $table               = $arg_ref->{table};
  my $internal_id         = $arg_ref->{internal_id};

  my $dep_sth         = $self->get_dep_sth($dbi);
  my $ins_go_dep_sth  = $self->get_dep_go_sth($dbi);
  my $ins_ox_sth2     = $self->get_add_dep_ox($dbi);
  my $ins_ix_sth      = $self->get_ins_ix_sth($dbi);
  my $get_object_xref_id_sth = $self->get_ox_id_master_sth($dbi);


  local $ins_ox_sth2->{RaiseError} = 0;  # want to see duplicates and not die automatically
  local $ins_ox_sth2->{PrintError} = 0;

  my $object_xref_id;
  while(my $master_xref_id = pop(@$master_xref_ids)){
    my ($dep_xref_id, $link);
    $dep_sth->execute($master_xref_id);
    $dep_sth->bind_columns(\$dep_xref_id, \$link);
    while($dep_sth->fetch){
      $ins_ox_sth2->execute($internal_id, $dep_xref_id, $table, 'DEPENDENT', $master_xref_id);
      $get_object_xref_id_sth->execute($internal_id, $dep_xref_id, $table, 'DEPENDENT', $master_xref_id);
      $object_xref_id = ($get_object_xref_id_sth->fetchrow_array())[0];
      if($ins_ox_sth2->err){
	my $err = $ins_ox_sth2->errstr;
	if($err =~ /Duplicate/){
	  $$duplicate_dep_count++;
	  next;
	}
	else{
	  die "Problem loading error is $err\n";
	} 
      }
      $ins_ix_sth->execute($object_xref_id);
      push @$master_xref_ids, $dep_xref_id; # get the dependent, dependents just in case

      if(defined($link) and $link ne ""){ # we have a go term linkage type
	$ins_go_dep_sth->execute($$object_xref_id, $link, $master_xref_id);
      }
    }
  }
  return;
}

1;

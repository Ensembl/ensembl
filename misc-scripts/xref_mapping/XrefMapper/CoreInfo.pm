package XrefMapper::CoreInfo;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use strict;
use warnings;
use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;

# Get info from the core database.

# Need to load tables:-
#
# gene_transcript_translation 
# gene_stable_id
# transcript_stable_id
# translation_stable_id





# Also for Direct xref process these and add to object_xref
#   also add dependents while we are at it.



sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->verbose($mapper->verbose);
  return $self;
}



sub get_core_data {
 my $self = shift;

# gene_transcript_translation 
# gene_stable_id
# transcript_stable_id
# translation_stable_id

  my $object_xref_id;
  my $ox_sth = $self->xref->dbc->prepare("select max(object_xref_id) from object_xref");
  $ox_sth->execute();
  $ox_sth->bind_columns(\$object_xref_id);
  $ox_sth->fetch();
  $ox_sth->finish;

 # load table gene_transcript_translation 

 my $ins_sth =  $self->xref->dbc->prepare("insert into gene_transcript_translation (gene_id, transcript_id, translation_id) values (?, ?, ?)"); 

 my $sql = "select tn.gene_id, tn.transcript_id, tl.translation_id from transcript tn left join translation tl on tl.transcript_id = tn.transcript_id";
 my $sth = $self->core->dbc->prepare($sql);
 $sth->execute();
 my  ($gene_id, $transcript_id, $translation_id);
 $sth->bind_columns(\$gene_id, \$transcript_id, \$translation_id); 
 while($sth->fetch()){
   $ins_sth->execute($gene_id, $transcript_id, $translation_id);
 }
 $ins_sth->finish;
 $sth->finish;


 # load table xxx_stable_id
 my ($id, $stable_id);
 foreach my $table (qw(gene transcript translation)){
   my $sth = $self->core->dbc->prepare("select ".$table."_id, stable_id from ".$table."_stable_id");
   my $ins_sth = $self->xref->dbc->prepare("insert into ".$table."_stable_id (internal_id, stable_id) values(?, ?)");
   $sth->execute();
   $sth->bind_columns(\$id, \$stable_id);
   while($sth->fetch){
     $ins_sth->execute($id, $stable_id);
   }
   $ins_sth->finish;
   $sth->finish;
 }

 $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('core_data_loaded',now())");
 $sth->execute();
 $sth->finish;

 # Now process the direct xrefs and add data to the object xrefs remember dependent xrefs.

 my $ins_ox_sth = $self->xref->dbc->prepare("insert into object_xref (object_xref_id, ensembl_id, xref_id, ensembl_object_type, linkage_type) values(?, ?, ?, ?, ?)");


 # Direct xrefs can be considered to be 100% matching
 my $ins_ix_sth = $self->xref->dbc->prepare("insert into identity_xref (object_xref_id, query_identity, target_identity) values(?, 100, 100)");

 local $ins_ox_sth->{RaiseError};  # want to see duplicates and not add de

 local $ins_ox_sth->{PrintError}; 
 

 my $ins_go_sth = $self->xref->dbc->prepare("insert into go_xref (object_xref_id, linkage_type, source_xref_id) values(?,?,?)");
 my $dep_sth    = $self->xref->dbc->prepare("select dependent_xref_id, linkage_annotation from dependent_xref where master_xref_id = ?");

my $stable_sql=(<<SQL);
  SELECT so.name, dx.general_xref_id, s.internal_id, dx.ensembl_stable_id 
    FROM source so, xref x, TYPE_direct_xref dx left join TYPE_stable_id s on s.stable_id = dx.ensembl_stable_id
      WHERE x.xref_id = dx.general_xref_id and x.source_id = so.source_id 
SQL
                      
 my %err_count;
 
 foreach my $table (qw(gene transcript translation)){
   my ($dbname, $xref_id, $internal_id, $stable_id);
   my $sql = $stable_sql;
   $sql =~ s/TYPE/$table/g;
   my $sth = $self->xref->dbc->prepare($sql);
#   print "sql = $sql\n";
   $sth->execute();
   $sth->bind_columns(\$dbname, \$xref_id, \$internal_id, \$stable_id);
   my $count =0;
   my $duplicate_direct_count = 0;
   my $duplicate_dependent_count = 0;
   while($sth->fetch){
     if(!defined($internal_id)){ # not found either it is an internal id already or stable_id no longer exists
       if($stable_id =~ /^\d+$/){
          $internal_id = $stable_id;
       }
       else{
	 if(!defined($err_count{$dbname}) or $err_count{$dbname} < 10){
	   print "Could not find stable id $stable_id in table to get the internal id hence ignoring!!! (for $dbname)\n" if($self->verbose);
	 }
	 $err_count{$dbname}++;
#	 $err_count++;
         next;
       }
     }
     $object_xref_id++;
     $count++;
     my @master_xref_ids;
     $ins_ox_sth->execute($object_xref_id, $internal_id, $xref_id, $table, 'DIRECT');
     if($ins_ox_sth->err){
       $duplicate_direct_count++;
       next; #duplicate
     }
     else{
       $ins_ix_sth->execute($object_xref_id);
       push  @master_xref_ids, $xref_id;
     }
     while(my $master_xref_id = pop(@master_xref_ids)){
       my ($dep_xref_id, $link);
       $dep_sth->execute($master_xref_id);
       $dep_sth->bind_columns(\$dep_xref_id, \$link);
       while($dep_sth->fetch){
         $object_xref_id++;
         $ins_ox_sth->execute($object_xref_id, $internal_id, $dep_xref_id, $table, 'DEPENDENT');
	 if($ins_ox_sth->err){
	   my $err = $ins_ox_sth->errstr;
	   if($err =~ /Duplicate/){
	     $duplicate_dependent_count++;
	     next;
	   }
	   else{
	     die "Problem loading error is $err\n";
	   } 
	 }
	 $ins_ix_sth->execute($object_xref_id);
	 push @master_xref_ids, $dep_xref_id; # get the dependent, dependents just in case

         if(defined($link) and $link ne ""){ # we have a go term linkage type
           $ins_go_sth->execute($object_xref_id, $link, $master_xref_id);
         }
       }
     }
   }
   $sth->finish;
   if($duplicate_direct_count or $duplicate_dependent_count){
     print "duplicate entrys ignored for $duplicate_direct_count direct xrefs and  $duplicate_dependent_count dependent xrefs\n" if($self->verbose);
   }
#   if($err_count or $self->verbose){
#     print STDERR $count." direct_xrefs added to ensembl ".$table."s BUT $err_count stable ids could not be found\n";
#   }
 }
 foreach my $key (%err_count){
   print STDERR "*WARNING*: ".$err_count{$key}." direct xrefs for database $key could not be added as their stable_ids could not be found\n";
 }
 $ins_go_sth->finish;
 $ins_ox_sth->finish;
 $dep_sth->finish;

 $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('direct_xrefs_parsed',now())");
 $sth->execute();
 $sth->finish;



}

1;

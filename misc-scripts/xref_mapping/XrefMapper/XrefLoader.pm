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

package XrefMapper::XrefLoader;
use strict;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use warnings;
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
  $self->mapper($mapper);
  return $self;
}


sub mapper{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_mapper} = $arg );
  return $self->{_mapper};
}

sub update{
  my ($self, $arg) = @_;
  # remove xref, object_xref, identity_xref, depenedent_xref, go_xref, unmapped_object, (interpro???), external_synonym, projections.


 my $verbose = $self->mapper->verbose;
  my $core_dbi = $self->core->dbc;
  my $xref_dbi = $self->xref->dbc;

  #####################################
  # first remove all the projections. #
  #####################################
  print "Deleting all PROJECTIONs from this database\n" if $verbose;

  my $sql = "DELETE es FROM xref x, external_synonym es WHERE x.xref_id = es.xref_id and x.info_type = 'PROJECTION'";
  my $sth = $core_dbi->prepare($sql);
  my $affected_rows = $sth->execute();
  print "\tDeleted $affected_rows PROJECTED external_synonym row(s)\n" if $verbose;

  # Delete all ontologies, as they are done by a separate pipeline
  $sql = <<SQL;
DELETE ontology_xref, object_xref, xref, dependent_xref
FROM ontology_xref, object_xref, xref 
LEFT JOIN dependent_xref on xref_id = dependent_xref_id
WHERE ontology_xref.object_xref_id = object_xref.object_xref_id AND object_xref.xref_id = xref.xref_id
SQL
  $sth = $core_dbi->prepare($sql);
  $affected_rows = $sth->execute();
  print "\tDeleted $affected_rows PROJECTED ontology_xref row(s)\n" if $verbose;

  $sql = "DELETE object_xref FROM object_xref, xref WHERE object_xref.xref_id = xref.xref_id AND xref.info_type = 'PROJECTION'";
  $sth = $core_dbi->prepare($sql);
  $affected_rows = $sth->execute();
  print "\tDeleted $affected_rows PROJECTED object_xref row(s)\n" if $verbose;

  $sql = "DELETE xref FROM xref WHERE xref.info_type = 'PROJECTION'";
  $sth = $core_dbi->prepare($sql);
  $affected_rows = $sth->execute();
  print "\tDeleted $affected_rows PROJECTED xref row(s)\n" if $verbose;

  $sth->finish;

  #########################################
  # Get source_id to external_db_id       #
  #########################################

  my %name_to_external_db_id;
  $sql = "select external_db_id, db_name from external_db";
  $sth = $core_dbi->prepare($sql);
  $sth->execute();
  my ($id, $name);
  $sth->bind_columns(\$id, \$name);
  while($sth->fetch()){
    $name_to_external_db_id{$name} = $id;
   }
  $sth->finish;

  my %source_id_to_external_db_id;
 
  $sql = 'select s.source_id, s.name from source s, xref x where x.source_id = s.source_id group by s.source_id'; # only get those of interest
  $sth = $xref_dbi->prepare($sql);
  $sth->execute();
  $sth->bind_columns(\$id, \$name);
  while($sth->fetch()){
    if(defined($name_to_external_db_id{$name})){
      $source_id_to_external_db_id{$id} = $name_to_external_db_id{$name};
    }
    elsif( $name =~ /notransfer$/){
    }	
    else{
      die "ERROR: Could not find $name in external_db table please add this too continue";
    }
  }
  $sth->finish;


  $sth = $xref_dbi->prepare("update xref set dumped = null where dumped != 'NO_DUMP_ANOTHER_PRIORITY'"); # just incase this is being ran again
  $sth->execute;
  $sth->finish;



  
  ######################################
  # For each external_db to be updated #
  # Delete the existing ones           # 
  ######################################
  my ($count);
  $sth = $xref_dbi->prepare('select s.name, count(*) from xref x, object_xref ox, source s where ox.xref_id = x.xref_id and x.source_id = s.source_id group by s.name');
  $sth->execute();
  $sth->bind_columns(\$name,\$count);

  my $synonym_sth  =  $core_dbi->prepare('DELETE external_synonym FROM external_synonym, xref WHERE external_synonym.xref_id = xref.xref_id AND xref.external_db_id = ?');
  my $go_sth       =  $core_dbi->prepare('DELETE ontology_xref.* FROM ontology_xref, object_xref, xref WHERE ontology_xref.object_xref_id = object_xref.object_xref_id AND object_xref.xref_id = xref.xref_id  AND xref.external_db_id = ?');
  my $identity_sth =  $core_dbi->prepare('DELETE identity_xref FROM identity_xref, object_xref, xref WHERE identity_xref.object_xref_id = object_xref.object_xref_id AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = ?');
  my $object_sth   =  $core_dbi->prepare('DELETE object_xref FROM object_xref, xref WHERE object_xref.xref_id = xref.xref_id AND xref.external_db_id = ?');
  my $dependent_sth = $core_dbi->prepare('DELETE d FROM dependent_xref d, xref x WHERE d.dependent_xref_id = x.xref_id and x.external_db_id = ?');
  my $xref_sth     =  $core_dbi->prepare('DELETE FROM xref WHERE xref.external_db_id = ?');
  my $unmapped_sth =  $core_dbi->prepare('DELETE FROM unmapped_object WHERE type="xref" and external_db_id = ?');


  my $transaction_start_sth  =  $core_dbi->prepare('start transaction');
  my $transaction_end_sth    =  $core_dbi->prepare('commit');

#
# ?? Is it faster to delete them all in one go with a external_db_id in (....) ???
# alternative load ottt etc that are not obtained from xrefs into xref table and then delete tables fully??
#



#  my $test =1;  # Can take a while so make optional when testing
#  if(!$test){
  $transaction_start_sth->execute();
  while($sth->fetch()){
    if(!defined($name_to_external_db_id{$name})){
      next;  #must end in notransfer
    }

    my $ex_id = $name_to_external_db_id{$name};

    print "Deleting data for $name from core before updating from new xref database\n" if ($verbose);
    $affected_rows = $synonym_sth->execute($ex_id);
    print "\tDeleted $affected_rows external_synonym row(s)\n" if $verbose;
    $affected_rows = $go_sth->execute($ex_id);
    print "\tDeleted $affected_rows ontology_xref row(s)\n" if $verbose;
    $affected_rows = $identity_sth->execute($ex_id);
    print "\tDeleted $affected_rows identity_xref row(s)\n" if $verbose;
    $affected_rows = $object_sth->execute($ex_id);  
    print "\tDeleted $affected_rows object_xref row(s)\n" if $verbose;
    $affected_rows = $dependent_sth->execute($ex_id);
    print "\tDeleted $affected_rows dependent_xref row(s)\n" if $verbose;
    $affected_rows = $xref_sth->execute($ex_id);
    print "\tDeleted $affected_rows xref row(s)\n" if $verbose;
    $affected_rows = $unmapped_sth->execute($ex_id);
    print "\tDeleted $affected_rows unmapped_object row(s)\n" if $verbose;
  }
  $sth->finish;
  $transaction_end_sth->execute();
#}
  $synonym_sth->finish;
  $go_sth->finish;  
  $identity_sth->finish;
  $object_sth->finish;  
  $dependent_sth->finish;
  $xref_sth->finish;
  $unmapped_sth->finish; 



  ##########################################
  # Get the offsets for object_xref, xref  #
  ##########################################

  $sth = $core_dbi->prepare('select MAX(xref_id) from xref');
  my $xref_offset;
  $sth->execute;
  $sth->bind_columns(\$xref_offset);
  $sth->fetch();
  $sth->finish;
  $xref_offset = 0 if(!defined($xref_offset));

  $self->add_meta_pair("xref_offset", $xref_offset);

  $sth = $core_dbi->prepare('select MAX(object_xref_id) from object_xref');
  my $object_xref_offset;
  $sth->execute;
  $sth->bind_columns(\$object_xref_offset);
  $sth->fetch();
  $sth->finish;
  $object_xref_offset = 0 if(!defined($object_xref_offset));

  $self->add_meta_pair("object_xref_offset", $object_xref_offset);


  ####################
  # Get analysis id's 
  ####################

  my %analysis_ids = $self->get_analysis();
  my $checksum_analysis_id; #do not populate until we know we need this


  print "xref offset is $xref_offset, object_xref offset is $object_xref_offset\n" if ($verbose);

  #####################################
  # Now add the new ones              #
  #####################################

     ###########################
     # SQL to get data from xref
     ###########################

  my $seq_sql =(<<DIRS);
SELECT x.xref_id, x.accession, x.label, x.version, x.description, x.info_text,
       ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type,
       i.query_identity, i.target_identity, i.hit_start, i.hit_end,
       i.translation_start, i.translation_end, i.cigar_line, i.score, i.evalue
  FROM xref x, object_xref ox, identity_xref i
    WHERE ox.ox_status = "DUMP_OUT" AND
          i.object_xref_id = ox.object_xref_id AND
          ox.xref_id = x.xref_id AND
          x.source_id = ? AND
          x.info_type = ? order by x.xref_id
DIRS

     my $seq_sth = $xref_dbi->prepare($seq_sql);


    ###########################
     # SQL to get data from xref without identity xref
     ###########################


  my $dir_sql =(<<DIRS);
SELECT x.xref_id, x.accession, x.label, x.version, x.description, x.info_text,
       ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type
  FROM xref x, object_xref ox
    WHERE ox.ox_status = "DUMP_OUT" AND
          ox.xref_id = x.xref_id AND
          x.source_id = ? AND
          x.info_type = ? order by x.xref_id
DIRS

     my $dir_sth = $xref_dbi->prepare($dir_sql);
 
#     $dependent_sth = $xref_dbi->prepare('select  x.xref_id, x.accession, x.label, x.version, x.description, x.info_text, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, d.master_xref_id from xref x, object_xref ox,  dependent_xref d where ox.ox_status = "DUMP_OUT" and ox.xref_id = x.xref_id and d.object_xref_id = ox.object_xref_id and x.source_id = ? and x.info_type = ? order by x.xref_id, ox.ensembl_id');
 
  my $dep_sql =(<<DSQL);
SELECT  x.xref_id, x.accession, x.label, x.version, x.description, x.info_text,
        ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, ox.master_xref_id 
   FROM xref x, object_xref ox 
     WHERE ox.ox_status = "DUMP_OUT" and 
           ox.xref_id = x.xref_id and 
           x.source_id = ? and 
           x.info_type = ? 
     ORDER BY x.xref_id, ox.ensembl_id
DSQL

 $dependent_sth = $xref_dbi->prepare($dep_sql);


  my $go_sql =(<<GSQL);
  SELECT  x.xref_id, x.accession, x.label, x.version, x.description, x.info_text, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, ox.master_xref_id, g.linkage_type,
       i.query_identity, i.target_identity, i.hit_start, i.hit_end,
       i.translation_start, i.translation_end, i.cigar_line, i.score, i.evalue
    FROM (xref x, object_xref ox, go_xref g, identity_xref i)
      WHERE ox.ox_status = "DUMP_OUT" and
            i.object_xref_id = ox.object_xref_id AND
            g.object_xref_id = ox.object_xref_id and
            x.xref_id = ox.xref_id and
            x.source_id = ? and x.info_type = ?
            order by x.xref_id, ox.ensembl_id
GSQL

     $go_sth = $xref_dbi->prepare($go_sql);

  my $go_count_sql = (<<GCNTSQL);
  SELECT  count(*)
    FROM (xref x, object_xref ox, go_xref g)
      WHERE ox.ox_status = "DUMP_OUT" and
            g.object_xref_id = ox.object_xref_id and
            x.xref_id = ox.xref_id and
            x.source_id = ? and x.info_type = ?
GCNTSQL
   
     # SQL to add data to core
     #########################
 
     my $add_identity_xref_sth  = $core_dbi->prepare('insert ignore into identity_xref (object_xref_id, xref_identity, ensembl_identity, xref_start, xref_end, ensembl_start, ensembl_end, cigar_line, score, evalue) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
     my $add_go_xref_sth        = $core_dbi->prepare('insert ignore into ontology_xref (object_xref_id, source_xref_id, linkage_type) values (?, ?, ?)');
     my $add_dependent_xref_sth = $core_dbi->prepare('insert ignore into dependent_xref (object_xref_id, master_xref_id, dependent_xref_id) values (?, ?, ?)');
     my $add_syn_sth            = $core_dbi->prepare('insert ignore into external_synonym (xref_id, synonym) values (?, ?)');
     my $add_release_info_sth   = $core_dbi->prepare('update external_db set db_release = ? where external_db_id = ?');

  $sth = $xref_dbi->prepare('select s.name, s.source_id, count(*), x.info_type, s.priority_description, s.source_release from xref x, object_xref ox, source s where ox.xref_id = x.xref_id  and x.source_id = s.source_id and ox_status = "DUMP_OUT" group by s.name, s.source_id, x.info_type');
  $sth->execute();
  my ($type, $source_id, $where_from, $release_info);
  $sth->bind_columns(\$name,\$source_id, \$count, \$type, \$where_from, \$release_info);
 
  $transaction_start_sth->execute();

  while($sth->fetch()){

    next if(!defined($name_to_external_db_id{$name}));

    if(defined($where_from) and $where_from ne ""){
      $where_from = "Generated via $where_from";
    }	
    my $ex_id = $name_to_external_db_id{$name};

    print "updating ($source_id) $name in core (for $type xrefs)\n" if ($verbose);

    my @xref_list=();  # process at end. Add synonyms and set dumped = 1;

    my $go_count_sth = $xref_dbi->prepare($go_count_sql);
    $go_count_sth->execute($source_id, $type);
    my ($go_data_present) = $go_count_sth->fetchrow_array;
    $go_count_sth->finish;
   
    # dump SEQUENCE_MATCH, DEPENDENT, DIRECT, COORDINATE_OVERLAP, INFERRED_PAIR, (MISC?? same as direct come from official naming)  

    ### If DIRECT ,         xref, object_xref,                  (order by xref_id)  # maybe linked to more than one?
    ### if INFERRED_PAIR    xref, object_xref
    ### if MISC             xref, object_xref 

    
    if($type eq "DIRECT" or $type eq "INFERRED_PAIR" or $type eq "MISC"){
      if ($go_data_present) {
       my $count = 0;
       $go_sth->execute($source_id, $type);
       my ($xref_id, $acc, $label, $version, $desc, $info,  $object_xref_id, $ensembl_id, $ensembl_type, $master_xref_id, $linkage_type); 
       my ( $query_identity, $target_identity, $hit_start, $hit_end, $translation_start, $translation_end, $cigar_line, $score, $evalue);
       $go_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id, \$ensembl_id, \$ensembl_type, \$master_xref_id, \$linkage_type,
                             \$query_identity, \$target_identity, \$hit_start, \$hit_end, \$translation_start, \$translation_end, \$cigar_line, \$score, \$evalue);
       my $last_xref = 0;
       while($go_sth->fetch){
	 if($last_xref != $xref_id){
	   push @xref_list, $xref_id;
	   $count++;
	   $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label, $version, $desc, $type, $info || $where_from, $core_dbi);
	   $last_xref = $xref_id;
	 }
         $object_xref_id = $self->add_object_xref($object_xref_offset, $object_xref_id, $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $analysis_ids{$ensembl_type}, $core_dbi);
         $add_go_xref_sth->execute( ($object_xref_id+$object_xref_offset), 0, $linkage_type);
         $add_identity_xref_sth->execute( ($object_xref_id+$object_xref_offset), $query_identity, $target_identity, $hit_start, $hit_end,
                                         $translation_start, $translation_end, $cigar_line, $score, $evalue) if $translation_start;
       }
       print "Direct GO $count\n" if ($verbose);
     }
     else{
       my $count = 0;
       $seq_sth->execute($source_id, $type);
       my ($xref_id, $acc, $label, $version, $desc, $info, $object_xref_id, $ensembl_id, $ensembl_type); 
       my ( $query_identity, $target_identity, $hit_start, $hit_end, $translation_start, $translation_end, $cigar_line, $score, $evalue);
       $seq_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id, \$ensembl_id, \$ensembl_type,
                             \$query_identity, \$target_identity, \$hit_start, \$hit_end, \$translation_start, \$translation_end, \$cigar_line, \$score, \$evalue);
       my $last_xref = 0;
       while($seq_sth->fetch){
         if($last_xref != $xref_id){
	  push @xref_list, $xref_id;
	  $count++;
	  $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label, $version, $desc, $type, $info || $where_from, $core_dbi);
	  $last_xref = $xref_id;
        }
        $object_xref_id = $self->add_object_xref($object_xref_offset, $object_xref_id, $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $analysis_ids{$ensembl_type}, $core_dbi);
        $add_identity_xref_sth->execute( ($object_xref_id+$object_xref_offset), $query_identity, $target_identity, $hit_start, $hit_end,
                                         $translation_start, $translation_end, $cigar_line, $score, $evalue) if $translation_start;
      }  
      print "DIRECT $count\n" if ($verbose);
     }
    }
    ### IF CHECKSUM,        xref, object_xref
    # 1:m mapping between object & xref
    elsif($type eq 'CHECKSUM') {
      #If we had a checksum then get the analysis. Avoids unecessary analysis entries
      if(! defined $checksum_analysis_id) {
        $checksum_analysis_id = $self->get_single_analysis('xrefchecksum');
      }
      my $count = 0;
      $dir_sth->execute($source_id, $type);
      my $last_xref = 0;
      while(my $row = $dir_sth->fetchrow_arrayref()) {
        my ($xref_id, $acc, $label, $version, $desc, $info, $object_xref_id, $ensembl_id, $ensembl_type) = @{$row};
        if($last_xref != $xref_id) {
          push @xref_list, $xref_id;
          $count++;
          $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label, $version, $desc, $type, $info || $where_from, $core_dbi);
          $last_xref = $xref_id;
        }
        $object_xref_id = $self->add_object_xref($object_xref_offset, $object_xref_id, $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $checksum_analysis_id, $core_dbi);
      }
      print "CHECKSUM $count\n" if ($verbose);
    }
 
    ### If DEPENDENT,       xref, object_xref , dependent_xref  (order by xref_id)  # maybe linked to more than one?
 
   elsif($type eq "DEPENDENT"){
     if ($go_data_present) {
       my $count = 0;
       $go_sth->execute($source_id, $type);
       my ($xref_id, $acc, $label, $version, $desc, $info,  $object_xref_id, $ensembl_id, $ensembl_type, $master_xref_id, $linkage_type);
       my ( $query_identity, $target_identity, $hit_start, $hit_end, $translation_start, $translation_end, $cigar_line, $score, $evalue);
       $go_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id, \$ensembl_id, \$ensembl_type, \$master_xref_id, \$linkage_type,
                             \$query_identity, \$target_identity, \$hit_start, \$hit_end, \$translation_start, \$translation_end, \$cigar_line, \$score, \$evalue);
       my $last_xref = 0;
       while($go_sth->fetch){
	 if($last_xref != $xref_id){
	   push @xref_list, $xref_id;
	   $count++;
	   $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label, $version, $desc, $type, $info || $where_from, $core_dbi);
	   $last_xref = $xref_id;
	 }
         $object_xref_id = $self->add_object_xref($object_xref_offset, $object_xref_id, $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $analysis_ids{$ensembl_type}, $core_dbi);
	 if(defined($master_xref_id)){  # need to sort this out as all should habe one really. (interpro generates go without these!!)
	   $add_dependent_xref_sth->execute(($object_xref_id+$object_xref_offset), ($master_xref_id+$xref_offset), ($xref_id+$xref_offset) );
	   $add_go_xref_sth->execute( ($object_xref_id+$object_xref_offset), ($master_xref_id+$xref_offset), $linkage_type);
           $add_identity_xref_sth->execute( ($object_xref_id+$object_xref_offset), $query_identity, $target_identity, $hit_start, $hit_end,
                                         $translation_start, $translation_end, $cigar_line, $score, $evalue) if $translation_start;
	 }
	 else {
	     $add_go_xref_sth->execute( ($object_xref_id+$object_xref_offset), 0, $linkage_type);
	 }
       }       
       print "GO $count\n" if ($verbose);     
     }
     else{
       my $count = 0;
       my $ox_count = 0;
       my @master_problems;
       my $err_master_count=0;
       $dependent_sth->execute($source_id, $type);
       my ($xref_id, $acc, $label, $version, $desc, $info, $object_xref_id, $ensembl_id, $ensembl_type, $master_xref_id); 
       $dependent_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id, \$ensembl_id, \$ensembl_type, \$master_xref_id);
       my $last_xref = 0;
       my $last_ensembl = 0;
       while($dependent_sth->fetch){
	 if($last_xref != $xref_id){
	   push @xref_list, $xref_id;
	   $count++;
	   $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label || $acc, $version, $desc, $type, $info || $where_from, $core_dbi);
	 }
	 if($last_xref != $xref_id or $last_ensembl != $ensembl_id){
	   $object_xref_id = $self->add_object_xref($object_xref_offset, $object_xref_id, $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $analysis_ids{$ensembl_type}, $core_dbi);
	   if (defined($master_xref_id)){ # need to sort this out for FlyBase since there are EMBL direct entries from the GFF and dependent xrefs from Uniprot
	     $add_dependent_xref_sth->execute(($object_xref_id+$object_xref_offset), ($master_xref_id+$xref_offset), ($xref_id+$xref_offset) );
	   }
	   else{
	     if($err_master_count < 10){
	       push @master_problems, $acc;
	     }
	     $err_master_count++;
	   }
	   $ox_count++;
	 }
	 $last_xref = $xref_id;
	 $last_ensembl = $ensembl_id;
       }
       if(@master_problems){
	 print "WARNING:: for $name $err_master_count problem master xrefs\nExamples are :-\t";
	 print join ", ",@master_problems;
	 print "\n";
       }
       print "DEP $count xrefs, $ox_count object_xrefs\n" if ($verbose);
     }
   }
   ### If SEQUENCE_MATCH   xref, object_xref,  identity_xref   (order by xref_id)  # maybe linked to more than one?

    elsif($type eq "SEQUENCE_MATCH"){
      my $count = 0;
      $seq_sth->execute($source_id, $type);
      my ($xref_id, $acc, $label, $version, $desc, $info, $object_xref_id, $ensembl_id, $ensembl_type); 
      my ( $query_identity, $target_identity, $hit_start, $hit_end, $translation_start, $translation_end, $cigar_line, $score, $evalue);
      $seq_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id, \$ensembl_id, \$ensembl_type,
			     \$query_identity, \$target_identity, \$hit_start, \$hit_end, \$translation_start, \$translation_end, \$cigar_line, \$score, \$evalue);
      my $last_xref = 0;
      while($seq_sth->fetch){
        if($last_xref != $xref_id){
	  push @xref_list, $xref_id;
	  $count++;
	  $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label, $version, $desc, $type, $info || $where_from, $core_dbi);
	  $last_xref = $xref_id;
        }
        $object_xref_id = $self->add_object_xref ($object_xref_offset, $object_xref_id, $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $analysis_ids{$ensembl_type}, $core_dbi);
	$add_identity_xref_sth->execute( ($object_xref_id+$object_xref_offset), $query_identity, $target_identity, $hit_start, $hit_end, 
					 $translation_start, $translation_end, $cigar_line, $score, $evalue);  
      }  
      print "SEQ $count\n" if ($verbose);
    }
    else{
      print "PROBLEM:: what type is $type\n";
    }	


    # Transfer data for synonym and set xref database xrefs to dumped.
    if(@xref_list){
      my $syn_count = 0;
      my $syn_sql = "select xref_id, synonym from synonym where xref_id in(".join(", ",@xref_list).")";
      my $syn_sth    = $xref_dbi->prepare($syn_sql);
      $syn_sth->execute();
    
      my ($xref_id, $syn);
      $syn_sth->bind_columns(\$xref_id, \$syn);
      while($syn_sth->fetch()){
	$add_syn_sth->execute(($xref_id+$xref_offset), $syn);
	$syn_count++;
      }
      $syn_sth->finish;

      print "\tadded $syn_count synonyms\n" if($syn_count);
      my $xref_dumped_sth = $xref_dbi->prepare("update xref set dumped = 'MAPPED' where xref_id in (".join(", ",@xref_list).")");
      $xref_dumped_sth->execute() || die "Could not set dumped status"; 
      $xref_dumped_sth->finish;
    }	
 
    # Update the core databases release in for source form the xref database
    if(defined($release_info) and $release_info ne "1"){
       $add_release_info_sth->execute($release_info, $ex_id) || die "Failed to add release info **$release_info** for external source $ex_id\n";
    }
  }
  $sth->finish;
  $seq_sth->finish;
  $dir_sth->finish;
  $transaction_end_sth->execute();


  #######################################
  # Remember to do unmapped entries
  # 1) make sure the reason exist/create them and get the ids for these.
  # 2) Process where dumped is null and type = DIRECT, DEPENDENT, SEQUENCE_MATCH, MISC seperately
  ########################################
  my %summary_failed;
  my %desc_failed;
  my %reason_id;

  # Get the cutoff values
  $sth = $xref_dbi->prepare("select distinct s.name, m.percent_query_cutoff, m.percent_target_cutoff from source s, source_mapping_method sm, mapping m where sm.source_id = s.source_id and sm.method = m.method");
  $sth->execute();
  my ($source_name, $q_cut, $t_cut);
  $sth->bind_columns(\$source_name, \$q_cut, \$t_cut);

  while ($sth->fetch) {
    $summary_failed{$source_name} = "Failed to match at thresholds";
    $desc_failed{$source_name}    = "Unable to match at the thresholds of $q_cut\% for the query or $t_cut\% for the target";
  }
  $sth->finish;

  $summary_failed{"NO_STABLE_ID"} = "Failed to find Stable ID";
  $desc_failed{"NO_STABLE_ID"}    = "Stable ID that this xref was linked to no longer exists";

  $summary_failed{"FAILED_MAP"} = "Failed to match";
  $desc_failed{"FAILED_MAP"}    = "Unable to match to any ensembl entity at all";

  $summary_failed{"NO_MAPPING"} = "No mapping done";
  $desc_failed{"NO_MAPPING"}    = "No mapping done for this type of xref";

  $summary_failed{"MASTER_FAILED"} = "Master failed";
  $desc_failed{"MASTER_FAILED"}    = "The dependent xref was not matched due to the master xref not being mapped";

  $summary_failed{"NO_MASTER"} = "No Master";
  $desc_failed{"NO_MASTER"}    = "The dependent xref was not matched due to there being no master xref";

  

  foreach my $key (keys %desc_failed){
    $sth = $core_dbi->prepare("select unmapped_reason_id from unmapped_reason where full_description like '".$desc_failed{$key}."'");
    $sth->execute();
    my $failed_id=undef;
    $sth->bind_columns(\$failed_id);
    $sth->fetch;
    $sth->finish;
    if(!defined($failed_id)){
      $sth = $core_dbi->prepare('insert into unmapped_reason (summary_description, full_description) values("'.$summary_failed{$key}.'", "'.$desc_failed{$key}.'")');
      $sth->execute();
      $failed_id = $sth->{'mysql_insertid'};
      $sth->finish
    }
    $reason_id{$key} = $failed_id;
  }

  $transaction_start_sth->execute();

  ##########
  # DIRECT #
  ##########

  my $dbname;
  $sql =(<<DIR);
  SELECT  x.xref_id, x.accession, x.version, x.label, x.description, x.info_type, x.info_text, s.name 
    FROM source s,xref x
      LEFT JOIN  object_xref ox ON ox.xref_id = x.xref_id
      WHERE x.source_id = s.source_id 
        AND x.dumped is null 
        AND ox.ox_status != 'FAILED_PRIORITY'
        AND x.info_type = 'DIRECT'
DIR

  my $direct_unmapped_sth = $xref_dbi->prepare($sql);
  my ($xref_id, $acc, $version, $label, $desc, $info);
  $direct_unmapped_sth->execute();
  $direct_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname);

  my $set_unmapped_sth       =  $core_dbi->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id ) values ('xref', ?, ?, ?, ?)");

  my @xref_list = ();
  my $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
  while($direct_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    if(defined($name_to_external_db_id{$dbname})){
      $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label, $version, $desc, 'UNMAPPED', $info, $core_dbi);   
      $set_unmapped_sth->execute($analysis_id, $ex_id, $acc, $reason_id{"NO_STABLE_ID"});
      push @xref_list, $xref_id;
    }
  }
  $direct_unmapped_sth->finish;
  $set_unmapped_sth->finish;


  if(@xref_list){
    my $xref_dumped_sth = $xref_dbi->prepare("update xref set dumped = 'UNMAPPED_NO_STABLE_ID' where xref_id in (".join(", ",@xref_list).")");
    $xref_dumped_sth->execute(); 
    $xref_dumped_sth->finish;
  }

  ########
  # MISC #
  ########

  $sql =(<<MIS);
  SELECT  x.xref_id, x.accession, x.version, x.label, x.description, x.info_type, x.info_text, s.name 
    FROM xref x, source s 
      WHERE x.source_id = s.source_id 
        AND x.dumped is null 
        AND x.info_type = 'MISC'
MIS

  my $misc_unmapped_sth = $xref_dbi->prepare($sql);
  $misc_unmapped_sth->execute();
  $misc_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname);

  @xref_list = ();
  $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
  while($misc_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    if(defined($name_to_external_db_id{$dbname})){
      $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label, $version, $desc, 'UNMAPPED', $info, $core_dbi);   
      $set_unmapped_sth->execute($analysis_id, $ex_id, $acc, $reason_id{"NO_MAPPING"});
      push @xref_list, $xref_id;
    }	
  }
  $misc_unmapped_sth->finish;
  $set_unmapped_sth->finish;


  if(@xref_list){
    my $xref_dumped_sth = $xref_dbi->prepare("update xref set dumped = 'UNMAPPED_NO_MAPPING' where xref_id in (".join(", ",@xref_list).")");
    $xref_dumped_sth->execute(); 
    $xref_dumped_sth->finish;
  }

  #############
  # DEPENDENT #
  #############

  $sql = (<<DEP);
    SELECT  distinct x.xref_id, x.accession, x.version, x.label, x.description, x.info_type, x.info_text, s.name, mx.accession 
      FROM xref mx, source s, xref x 
          LEFT JOIN dependent_xref dx ON  dx.dependent_xref_id = x.xref_id
          LEFT JOIN object_xref ox ON ox.xref_id = x.xref_id
        WHERE x.source_id = s.source_id 
          AND dx.master_xref_id = mx.xref_id 
          AND x.dumped is null 
          AND ox.ox_status != 'FAILED_PRIORITY'
          AND x.info_type = 'DEPENDENT'
          ORDER BY s.name, x.accession
DEP

  my $dep_unmapped_sth = $xref_dbi->prepare($sql);
  $dep_unmapped_sth->execute();
  my $parent;
  $dep_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname, \$parent);

  $set_unmapped_sth  =  $core_dbi->prepare("insert ignore into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id, parent ) values ('xref', ?, ?, ?, '".$reason_id{"MASTER_FAILED"}."', ?)");

  @xref_list = ();
  my $last_acc= 0;
  while($dep_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    if(!defined($ex_id)){
      next;
    }
    if($last_acc ne $acc){
      $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label||$acc, $version, $desc, 'UNMAPPED', $info, $core_dbi);
    }
    $last_acc = $acc;
    $set_unmapped_sth->execute($analysis_id, $ex_id, $acc, $parent);
    push @xref_list, $xref_id;
  }
  $dep_unmapped_sth->finish;
  $set_unmapped_sth->finish;


  if(@xref_list){
    my $xref_dumped_sth = $xref_dbi->prepare("update xref set dumped = 'UNMAPPED_MASTER_FAILED' where xref_id in (".join(", ",@xref_list).")");
    $xref_dumped_sth->execute(); 
    $xref_dumped_sth->finish;
  }

  ##################
  # SEQUENCE_MATCH #
  ##################

  $sql = (<<SEQ);
    SELECT  x.xref_id, x.accession, x.version, x.label, x.description, x.info_type, x.info_text, 
            s.name, px.sequence_type, 
            ox.ensembl_object_type, ox.ensembl_id,
            ix.query_identity, ix.target_identity, ox.ox_status
      FROM source s, primary_xref px, xref x
        LEFT JOIN object_xref ox ON ox.xref_id = x.xref_id
        LEFT JOIN identity_xref ix ON ix.object_xref_id = ox.object_xref_id
      WHERE x.source_id = s.source_id
	  AND px.xref_id = x.xref_id
          AND x.dumped is null 
          AND x.info_type = 'SEQUENCE_MATCH'
          ORDER  BY x.xref_id
          
SEQ
# removed          AND ox.ox_status != 'FAILED_PRIORITY'

  my $seq_unmapped_sth = $xref_dbi->prepare($sql);
  $seq_unmapped_sth->execute();
  my ($ensembl_object_type, $ensembl_id, $q_id, $t_id, $seq_type, $status) ;
  $seq_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname, \$seq_type, \$ensembl_object_type, \$ensembl_id, \$q_id, \$t_id,\$status);

  my $set_unmapped_no_sth     = $core_dbi->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id, ensembl_object_type ) values ('xref', ?, ?, ?, '".$reason_id{"FAILED_MAP"}."', ?)");
  my $set_unmapped_failed_sth = $core_dbi->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id, query_score, target_score, ensembl_id, ensembl_object_type ) values ('xref', ?, ?, ?, ?,?,?,?,?)");


  @xref_list = ();
  my $last_xref = 0;
  my $unmapped_reason_id;
  while($seq_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    if(!defined($ex_id) or (defined($status) and $status eq "FAILED_PRIORITY") ){
      next;
    }
    if($last_xref != $xref_id){
      $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label, $version, $desc, 'UNMAPPED', $info, $core_dbi);
    }
    $last_xref = $xref_id;
    if(defined($ensembl_id)){
      $analysis_id= $analysis_ids{$ensembl_object_type};
      $unmapped_reason_id = $reason_id{$dbname};
      $set_unmapped_failed_sth->execute($analysis_id, $ex_id, $acc, $unmapped_reason_id, $q_id, $t_id, $ensembl_id, $ensembl_object_type );
    }
    else{
      if($seq_type eq "dna"){
	$ensembl_object_type = "Transcript";
      }
      else{
	$ensembl_object_type = "Translation";
      }	
      $analysis_id = $analysis_ids{$ensembl_object_type};
      $set_unmapped_no_sth->execute($analysis_id, $ex_id, $acc, $ensembl_object_type);
    }
    push @xref_list, $xref_id;
  }
  $seq_unmapped_sth->finish;
  $set_unmapped_no_sth->finish;
  $set_unmapped_failed_sth->finish;


  if(@xref_list){
    my $xref_dumped_sth = $xref_dbi->prepare("update xref set dumped = 'UNMAPPED_NO_MAPPING' where xref_id in (".join(", ",@xref_list).")");
    $xref_dumped_sth->execute(); 
    $xref_dumped_sth->finish;
  }

  ###########################
  # WEL (What ever is left).#
  ###########################
  
  # These are those defined as dependent but the master never existed and the xref and their descriptions etc are loaded first
  # with the dependencys added later so did not know they had no masters at time of loading.
  # (e.g. EntrezGene, WikiGene, MIN_GENE, MIM_MORBID)

 $sql = (<<WEL);
    SELECT  distinct x.xref_id, x.accession, x.version, x.label, x.description, x.info_type, x.info_text, s.name
      FROM source s, xref x 
        WHERE x.source_id = s.source_id 
          AND x.dumped is null 
          AND x.info_type = 'DEPENDENT'
WEL

  
  my $wel_unmapped_sth = $xref_dbi->prepare($sql);
  $wel_unmapped_sth->execute();
  $wel_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname);

  $set_unmapped_sth  =  $core_dbi->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id) values ('xref', ?, ?, ?, '".$reason_id{"NO_MASTER"}."')");

  $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
  @xref_list = ();
  while($wel_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    if(!defined($ex_id)){
      next;
    }
    $xref_id = $self->add_xref($xref_offset, $xref_id, $ex_id, $acc, $label, $version, $desc, 'UNMAPPED', $info, $core_dbi);
    $set_unmapped_sth->execute($analysis_id, $ex_id, $acc);
    push @xref_list, $xref_id;
  }
  $wel_unmapped_sth->finish;
  $set_unmapped_sth->finish;

  if(@xref_list){
    my $xref_dumped_sth = $xref_dbi->prepare("update xref set dumped = 'UNMAPPED_NO_MASTER' where xref_id in (".join(", ",@xref_list).")");
    $xref_dumped_sth->execute(); 
    $xref_dumped_sth->finish;
  }


  $transaction_end_sth->execute();

  my $sth_stat = $xref_dbi->prepare("insert into process_status (status, date) values('core_loaded',now())");
  $sth_stat->execute();
  $sth_stat->finish;



}


sub get_analysis{
  my $self = shift;
  my %typeToLogicName = ( 'Gene'        => 'xrefexoneratedna',
                          'Transcript'  => 'xrefexoneratedna',
                          'Translation' => 'xrefexonerateprotein');
  my %analysis_id;
  foreach my $key (qw(Gene Transcript Translation)){
    my $logic_name = $typeToLogicName{$key};
    $analysis_id{$key} = $self->get_single_analysis($logic_name);    
  }
  return %analysis_id;
}

sub get_single_analysis {
  my ($self, $logic_name) = @_;
  my $h = $self->core->dbc()->sql_helper();
  my $analysis_ids = $h->execute_simple(
    -SQL => 'SELECT analysis_id FROM analysis WHERE logic_name=?', 
    -PARAMS => [$logic_name] 
  );
  my $analysis_id;
  
  if(@{$analysis_ids}) {
    $analysis_id = $analysis_ids->[0];
  }
  else {
    print "No analysis with logic_name $logic_name found, creating ...\n" if ($self->verbose);
    # TODO - other fields in analysis table
    $self->core()->dbc()->sql_helper()->execute_update(
      -SQL => 'INSERT INTO analysis (logic_name, created) VALUES (?,NOW())',
      -PARAMS => [$logic_name],
      -CALLBACK => sub {
        my ($sth) = @_;
        $analysis_id = $sth->{'mysql_insertid'};
        return;
      }
    );
  }
  
  return $analysis_id;
}


sub add_xref {
  my ($self, $offset, $xref_id, $external_db_id, $dbprimary_acc, $display_label, $version, $description, $info_type, $info_text, $dbc)  = @_;
  my $select_sth = $dbc->prepare("select xref_id from xref where dbprimary_acc = ? and external_db_id = ? and info_type = ? and info_text = ? and version = ?");
  my $insert_sth = $dbc->prepare("insert into xref (xref_id, external_db_id, dbprimary_acc, display_label, version, description, info_type, info_text) values (?, ?, ?, ?, ?, ?, ?, ?)");
  my $new_xref_id;
  $select_sth->execute($dbprimary_acc, $external_db_id, $info_type, $info_text, $version);
  $select_sth->bind_columns(\$new_xref_id);
  $select_sth->fetch();
  if (!$new_xref_id) {
    $insert_sth->execute(($xref_id+$offset), $external_db_id, $dbprimary_acc, $display_label, $version, $description, $info_type, $info_text);
    return $xref_id;
  } else {
    return $new_xref_id - $offset;
  }
}

sub add_object_xref {
  my ($self, $offset, $object_xref_id, $ensembl_id, $ensembl_object_type, $xref_id, $analysis_id, $dbc) = @_;
  my $select_sth = $dbc->prepare("select object_xref_id from object_xref where xref_id = ? and ensembl_object_type = ? and ensembl_id = ? and analysis_id = ?");
  my $insert_sth = $dbc->prepare("insert ignore into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id, analysis_id) values (?, ?, ?, ?, ?)");
  my $new_object_xref_id;
  $select_sth->execute($xref_id, $ensembl_object_type, $ensembl_id, $analysis_id);
  $select_sth->bind_columns(\$new_object_xref_id);
  $select_sth->fetch();
  if (!$new_object_xref_id) {
    $insert_sth->execute(($object_xref_id+$offset), $ensembl_id, $ensembl_object_type, $xref_id, $analysis_id);
    return $object_xref_id;
  } else {
    return $new_object_xref_id - $offset;
  }
}



1;

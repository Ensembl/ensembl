package XrefMapper::XrefLoader;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use strict;
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

  #####################################
  # first remove all the projections. #
  #####################################

  my $sql = "DELETE es FROM xref x, external_synonym es WHERE x.xref_id = es.xref_id and x.info_type = 'PROJECTION'";
  my $sth = $self->core->dbc->prepare($sql);
  $sth->execute();

  $sql = "DELETE object_xref FROM object_xref, xref WHERE object_xref.xref_id = xref.xref_id AND xref.info_type = 'PROJECTION'";
  $sth = $self->core->dbc->prepare($sql);
  $sth->execute();
  $sql = "DELETE xref FROM xref WHERE xref.info_type = 'PROJECTION'";
  $sth = $self->core->dbc->prepare($sql);
  $sth->execute();
  $sth->finish;

  #########################################
  # Get source_id to external_db_id       #
  #########################################

  my %name_to_external_db_id;
  $sql = "select external_db_id, db_name from external_db";
  $sth = $self->core->dbc->prepare($sql);
  $sth->execute();
  my ($id, $name);
  $sth->bind_columns(\$id, \$name);
  while($sth->fetch()){
    $name_to_external_db_id{$name} = $id;
   }
  $sth->finish;

  my %source_id_to_external_db_id;
 
  $sql = 'select s.source_id, s.name from source s, xref x where x.source_id = s.source_id group by s.source_id'; # only get those of interest
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();
  $sth->bind_columns(\$id, \$name);
  while($sth->fetch()){
     if(defined($name_to_external_db_id{$name})){
      $source_id_to_external_db_id{$id} = $name_to_external_db_id{$name};
    }
    else{
      die "ERROR: Could not find $name in external_db table please add this too continue\n";
    }
  }
  $sth->finish;


  # 100 is the dumped status for the dumpoed status for the "FAILED_PRIORIY"'s do not clear these.
  $sth = $self->xref->dbc->prepare("update xref set dumped = null where dumped != 100"); # just incase this is being ran again
  $sth->execute;
  $sth->finish;



  
  ######################################
  # For each external_db to be updated #
  # Delete the existing ones           # 
  ######################################

  $sth = $self->xref->dbc->prepare('select s.name, count(*) from xref x, object_xref ox, source s where ox.xref_id = x.xref_id and x.source_id = s.source_id and ox_status = "DUMP_OUT"  group by s.name');
  $sth->execute();
  my $count;
  $sth->bind_columns(\$name,\$count);

  my $synonym_sth  =  $self->core->dbc->prepare('DELETE external_synonym FROM external_synonym, xref WHERE external_synonym.xref_id = xref.xref_id AND xref.external_db_id = ?');
  my $go_sth       =  $self->core->dbc->prepare('DELETE FROM ontology_xref');
  my $identity_sth =  $self->core->dbc->prepare('DELETE identity_xref FROM identity_xref, object_xref, xref WHERE identity_xref.object_xref_id = object_xref.object_xref_id AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = ?');
  my $object_sth   =  $self->core->dbc->prepare('DELETE object_xref FROM object_xref, xref WHERE object_xref.xref_id = xref.xref_id AND xref.external_db_id = ?');
  my $dependent_sth = $self->core->dbc->prepare('DELETE d FROM dependent_xref d, xref x WHERE d.dependent_xref_id = x.xref_id and x.external_db_id = ?');
  my $xref_sth     =  $self->core->dbc->prepare('DELETE FROM xref WHERE xref.external_db_id = ?');
  my $unmapped_sth =  $self->core->dbc->prepare('DELETE FROM unmapped_object WHERE type="xref" and external_db_id = ?');


  my $transaction_start_sth  =  $self->core->dbc->prepare('start transaction');
  my $transaction_end_sth    =  $self->core->dbc->prepare('commit');

#
# ?? Is it faster to delete them all in one go with a external_db_id in (....) ???
# alternative load ottt etc that are not obtained from xrefs into xref table and then delete tables fully??
#



#  my $test =1;  # Can take a while so make optional when testing
#  if(!$test){
  $transaction_start_sth->execute();
  while($sth->fetch()){
    my $ex_id = $name_to_external_db_id{$name};

    print "Deleting data for $name from core before updating from new xref database\n" if ($verbose);
    $synonym_sth->execute($ex_id);
    if($name eq "GO"){
      $go_sth->execute();
    }
    $identity_sth->execute($ex_id);
    $object_sth->execute($ex_id);  
    $dependent_sth->execute($ex_id);
    $xref_sth->execute($ex_id);
    $unmapped_sth->execute($ex_id);
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

  $sth = $self->core->dbc->prepare('select MAX(xref_id) from xref');
  my $xref_offset;
  $sth->execute;
  $sth->bind_columns(\$xref_offset);
  $sth->fetch();
  $sth->finish;
  $xref_offset = 0 if(!defined($xref_offset));

  $self->add_meta_pair("xref_offset", $xref_offset);

  $sth = $self->core->dbc->prepare('select MAX(object_xref_id) from object_xref');
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

  my %analysis_ids = $self->get_analysis(); # 


  print "xref offset is $xref_offset, object_xref offset is $object_xref_offset\n" if ($verbose);

  #####################################
  # Now add the new ones              #
  #####################################

     ###########################
     # SQL to get data from xref
     ###########################

     my $direct_sth = $self->xref->dbc->prepare('select x.xref_id, x.accession, x.label, x.version, x.description, x.info_text, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type from xref x, object_xref ox  where ox.ox_status = "DUMP_OUT" and ox.xref_id = x.xref_id and x.source_id = ? and x.info_type = ? order by x.xref_id');
 
#     $dependent_sth = $self->xref->dbc->prepare('select  x.xref_id, x.accession, x.label, x.version, x.description, x.info_text, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, d.master_xref_id from xref x, object_xref ox,  dependent_xref d where ox.ox_status = "DUMP_OUT" and ox.xref_id = x.xref_id and d.object_xref_id = ox.object_xref_id and x.source_id = ? and x.info_type = ? order by x.xref_id, ox.ensembl_id');
 
    $dependent_sth = $self->xref->dbc->prepare('select  x.xref_id, x.accession, x.label, x.version, x.description, x.info_text, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, ox.master_xref_id from xref x, object_xref ox where ox.ox_status = "DUMP_OUT" and ox.xref_id = x.xref_id and x.source_id = ? and x.info_type = ? order by x.xref_id, ox.ensembl_id');


  my $go_sql =(<<GSQL);
  SELECT  x.xref_id, x.accession, x.label, x.version, x.description, x.info_text, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, ox.master_xref_id, g.linkage_type 
    FROM (xref x, object_xref ox, go_xref g) 
      WHERE ox.ox_status = "DUMP_OUT" and  
            g.object_xref_id = ox.object_xref_id and 
            x.xref_id = ox.xref_id and 
            x.source_id = ? and x.info_type = ? 
            order by x.xref_id, ox.ensembl_id
GSQL

     $go_sth = $self->xref->dbc->prepare($go_sql);

     my $seq_sth   =   $self->xref->dbc->prepare('select x.xref_id, x.accession, x.label, x.version, x.description, x.info_text, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, i.query_identity, i.target_identity, i.hit_start, i.hit_end, i.translation_start, i.translation_end, i.cigar_line, i.score, i.evalue from xref x, object_xref ox, identity_xref i  where ox.ox_status = "DUMP_OUT" and i.object_xref_id = ox.object_xref_id and ox.xref_id = x.xref_id and x.source_id = ? and x.info_type = ? order by x.xref_id');

     ########################
     # SQL to add data to core
     #########################
 
     my $add_xref_sth           = $self->core->dbc->prepare('insert into xref (xref_id, external_db_id, dbprimary_acc, display_label, version, description, info_type, info_text) values (?, ?, ?, ?, ?, ?, ?, ?)');
     my $add_object_xref_sth    = $self->core->dbc->prepare('insert into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id, analysis_id) values (?, ?, ?, ?, ?)');
     my $add_identity_xref_sth  = $self->core->dbc->prepare('insert into identity_xref (object_xref_id, xref_identity, ensembl_identity, xref_start, xref_end, ensembl_start, ensembl_end, cigar_line, score, evalue) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
     my $add_go_xref_sth        = $self->core->dbc->prepare('insert into ontology_xref (object_xref_id, linkage_type) values (?, ?)');
     my $add_dependent_xref_sth = $self->core->dbc->prepare('insert ignore into dependent_xref (object_xref_id, master_xref_id, dependent_xref_id) values (?, ?, ?)');
     my $add_syn_sth            = $self->core->dbc->prepare('insert ignore into external_synonym (xref_id, synonym) values (?, ?)');
     my $add_release_info_sth   = $self->core->dbc->prepare('update external_db set db_release = ? where external_db_id = ?');

  $sth = $self->xref->dbc->prepare('select s.name, s.source_id, count(*), x.info_type, s.priority_description, s.source_release from xref x, object_xref ox, source s where ox.xref_id = x.xref_id  and x.source_id = s.source_id and ox_status = "DUMP_OUT" group by s.name, s.source_id, x.info_type');
  $sth->execute();
  my ($type, $source_id, $where_from, $release_info);
  $sth->bind_columns(\$name,\$source_id, \$count, \$type, \$where_from, \$release_info);
 
  $transaction_start_sth->execute();

  while($sth->fetch()){
    if(defined($where_from) and $where_from ne ""){
      $where_from = "Generated via $where_from";
    }	
    my $ex_id = $name_to_external_db_id{$name};

    print "updating ($source_id) $name in core (for $type xrefs)\n" if ($verbose);

    my @xref_list=();  # process at end. Add synonyms and set dumped = 1;

   
    # dump SEQUENCE_MATCH, DEPENDENT, DIRECT, COORDINATE_OVERLAP, INFERRED_PAIR, (MISC?? same as direct come from official naming)  

    ### If DIRECT ,         xref, object_xref,                  (order by xref_id)  # maybe linked to more than one?
    ### if INFERRED_PAIR    xref, object_xref
    ### if MISC             xref, object_xref 

    
    if($type eq "DIRECT" or $type eq "INFERRED_PAIR" or $type eq "MISC"){
      my $count = 0;
      $direct_sth->execute($source_id, $type);
      my ($xref_id, $acc, $label, $version, $desc, $info, $object_xref_id, $ensembl_id, $ensembl_type); 
      $direct_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id, \$ensembl_id, \$ensembl_type);
      my $last_xref = 0;
      while($direct_sth->fetch){
        if($last_xref != $xref_id){
	  push @xref_list, $xref_id;
	  $count++;
	  $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type, $info || $where_from);
	  $last_xref = $xref_id;
        }
        $add_object_xref_sth->execute(($object_xref_id+$object_xref_offset), $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $analysis_ids{$ensembl_type});
      }  
      print "DIRECT $count\n" if ($verbose);
    }
 
    ### If DEPENDENT,       xref, object_xref , dependent_xref  (order by xref_id)  # maybe linked to more than one?
 
   elsif($type eq "DEPENDENT"){
     if($name eq "GO"){
       my $count = 0;
       $go_sth->execute($source_id, $type);
       my ($xref_id, $acc, $label, $version, $desc, $info,  $object_xref_id, $ensembl_id, $ensembl_type, $master_xref_id, $linkage_type); 
       $go_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id, \$ensembl_id, \$ensembl_type, \$master_xref_id, \$linkage_type);
       my $last_xref = 0;
       while($go_sth->fetch){
	 if($last_xref != $xref_id){
	   push @xref_list, $xref_id;
	   $count++;
	   $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type, $info || $where_from);
	   $last_xref = $xref_id;
	 }
	 if(defined($master_xref_id)){  # need to sort this out as all should habe one really. (interpro generates go without these!!)
	   $add_dependent_xref_sth->execute(($object_xref_id+$object_xref_offset), ($master_xref_id+$xref_offset), ($xref_id+$xref_offset) );
	 }
	 $add_object_xref_sth->execute(   ($object_xref_id+$object_xref_offset), $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $analysis_ids{$ensembl_type} );
	 $add_go_xref_sth->execute(       ($object_xref_id+$object_xref_offset), $linkage_type);
       }       
       print "GO $count\n" if ($verbose);     
     }
     else{
       my $count = 0;
       my $ox_count = 0;
       $dependent_sth->execute($source_id, $type);
       my ($xref_id, $acc, $label, $version, $desc, $info, $object_xref_id, $ensembl_id, $ensembl_type, $master_xref_id); 
       $dependent_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id, \$ensembl_id, \$ensembl_type, \$master_xref_id);
       my $last_xref = 0;
       my $last_ensembl = 0;
       while($dependent_sth->fetch){
	 if($last_xref != $xref_id){
	   push @xref_list, $xref_id;
	   $count++;
	   $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label || $acc, $version, $desc, $type, $info || $where_from);
	 }
	 if($last_xref != $xref_id or $last_ensembl != $ensembl_id){
	   $add_object_xref_sth->execute(($object_xref_id+$object_xref_offset), $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $analysis_ids{$ensembl_type});
	   $add_dependent_xref_sth->execute(($object_xref_id+$object_xref_offset), ($master_xref_id+$xref_offset), ($xref_id+$xref_offset) );
	   $ox_count++;
	 }
	 $last_xref = $xref_id;
	 $last_ensembl = $ensembl_id;
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
	  $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type, $info || $where_from);
	  $last_xref = $xref_id;
        }
        $add_object_xref_sth->execute(   ($object_xref_id+$object_xref_offset), $ensembl_id, $ensembl_type, ($xref_id+$xref_offset), $analysis_ids{$ensembl_type});
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
      my $syn_sql = "select xref_id, synonym from synonym where xref_id in(".join(", ",@xref_list).")";
      my $syn_sth    = $self->xref->dbc->prepare($syn_sql);
      $syn_sth->execute();
    
      my ($xref_id, $syn);
      $syn_sth->bind_columns(\$xref_id, \$syn);
      while($syn_sth->fetch()){
	$add_syn_sth->execute(($xref_id+$xref_offset), $syn)
      }
      $syn_sth->finish;

      my $xref_dumped_sth = $self->xref->dbc->prepare("update xref set dumped = 2 where xref_id in (".join(", ",@xref_list).")");
      $xref_dumped_sth->execute() || die "Could not set dumped status"; 
      $xref_dumped_sth->finish;
    }	
 
    # Update the core databases release in for source form the xref database
    if(defined($release_info) and $release_info ne "1"){
       $add_release_info_sth->execute($release_info, $ex_id) || die "Failed to add release info **$release_info** for external source $ex_id\n";
    }
  }
  $sth->finish;
  $transaction_end_sth->execute();


  #######################################
  # Remember to do unmapped entries
  # 1) make sure the reason exist/create them and get the ids for these.
  # 2) Process where dumped is null and type = DIRECT, DEPENDENT, SEQUENCE_MATCH, MISC seperately
  ########################################

  # Get the cutoff values
  $sth = $self->xref->dbc->prepare("select job_id, percent_query_cutoff, percent_target_cutoff from mapping limit 1");
  $sth->execute();
  my ($job_id, $q_cut, $t_cut);
  $sth->bind_columns(\$job_id, \$q_cut, \$t_cut);
  $sth->fetch;
  $sth->finish;


  my %summary_failed;
  my %desc_failed;
  my %reason_id;

  $summary_failed{"NO_STABLE_ID"} = "Failed to find Stable ID";
  $desc_failed{"NO_STABLE_ID"}    = "Stable ID that this xref was linked to no longer exists";

  $summary_failed{"FAILED_MAP"} = "Failed to match";
  $desc_failed{"FAILED_MAP"}    = "Unable to match to any ensembl entity at all";

  $summary_failed{"NO_MAPPING"} = "No mapping done";
  $desc_failed{"NO_MAPPING"}    = "No mapping done for this type of xref";

  $summary_failed{"FAILED_THRESHOLD"} = "Failed to match at thresholds";
  $desc_failed{"FAILED_THRESHOLD"}    = "Unable to match at the thresholds of $q_cut\% for the query or $t_cut\% for the target";

  $summary_failed{"MASTER_FAILED"} = "Master failed";
  $desc_failed{"MASTER_FAILED"}    = "The dependent xref was not matched due to the master xref not being mapped";

  $summary_failed{"NO_MASTER"} = "No Master";
  $desc_failed{"NO_MASTER"}    = "The dependent xref was not matched due to there being no master xref";

  

  foreach my $key (keys %desc_failed){
    $sth = $self->core->dbc->prepare("select unmapped_reason_id from unmapped_reason where full_description like '".$desc_failed{$key}."'");
    $sth->execute();
    my $failed_id=undef;
    $sth->bind_columns(\$failed_id);
    $sth->fetch;
    $sth->finish;
    if(!defined($failed_id)){
      $sth = $self->core->dbc->prepare('insert into unmapped_reason (summary_description, full_description) values("'.$summary_failed{$key}.'", "'.$desc_failed{$key}.'")');
      $sth->execute();
      $failed_id = $sth->{'mysql_insertid'};
      $sth->finish
    }
    $reason_id{$key} = $failed_id;
  }


  ## Dump interpro xrefs and interpro table
  # use NO_MAPPING as unmapped_reason
  
  # remove the old set
  # dump xrefs;
  # dump unmapped reasons
  # set xref status to dumped
  
  $transaction_start_sth->execute();

  #delete old set
  # 1 delete unmapped_object
  # 2 delete xrefs

  my $del_x_sth = $self->core->dbc->prepare('delete x from xref x, external_db e where x.external_db_id = e.external_db_id and e.db_name like "interpro"') || die "Could not prepare interpro xref deletion";
  $del_x_sth->execute() || die "Problem executing deletion of interpro xrefs";

  my $del_uo_sth = $self->core->dbc->prepare('delete x from unmapped_object x, external_db e where x.external_db_id = e.external_db_id and e.db_name like "interpro"') || die "Could not prepare interpro unmapped_object deletion";
  $del_uo_sth->execute() || die "Problem executing deletion of interpro xrefs";

  my $get_xref_interpro_sth  = $self->xref->dbc->prepare("select x.xref_id, x.accession, x.version, x.label, x.description, x.info_type, x.info_text from xref x ,source s where s.source_id = x.source_id and s.name like 'Interpro'");

  my $get_interpro_sth       = $self->xref->dbc->prepare("select interpro, pfam from interpro");
  my $add_interpro_sth       = $self->core->dbc->prepare("insert into interpro (interpro_ac, id) values (?, ?)");
  my $set_unmapped_sth       =  $self->core->dbc->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id ) values ('xref', ?, ?, ?, ?)");
  my @xref_list =();
  
  
  my $ex_id = $name_to_external_db_id{"Interpro"};
  my $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
 
  
  $get_xref_interpro_sth->execute();
  my ($xref_id, $acc, $version, $label, $desc, $info);
  $get_xref_interpro_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info);
  
  while($get_xref_interpro_sth->fetch){
    $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type, $info);
    $set_unmapped_sth->execute($analysis_id, $ex_id, $acc, $reason_id{"NO_MAPPING"} );
    push @xref_list, $xref_id;
  }
  $get_xref_interpro_sth->finish;
  $set_unmapped_sth->finish;

  if(@xref_list){
    my $xref_dumped_sth = $self->xref->dbc->prepare("update xref set dumped = 3 where xref_id in (".join(", ",@xref_list).")");
    $xref_dumped_sth->execute(); 
    $xref_dumped_sth->finish;
  }

  # delete all entries in interpro table
  my $del_sth = $self->core->dbc->prepare("delete from interpro");
  $del_sth->execute;
  $del_sth->finish;
  
  # add new entries to interpro table
  $get_interpro_sth->execute();
  my ($inter);
  $get_interpro_sth->bind_columns(\$inter,\$id);
  while($get_interpro_sth->fetch){
    $add_interpro_sth->execute($inter, $id);
  }
  $get_interpro_sth->finish;
  $add_interpro_sth->finish;

  $transaction_end_sth->execute();


#  foreach my $type (qw(MISC DEPENDENT DIRECT SEQUENCE_MATCH INFERRED_PAIR)){
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

  my $direct_unmapped_sth = $self->xref->dbc->prepare($sql);
  $direct_unmapped_sth->execute();
  $direct_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname);

  @xref_list = ();
  $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
  while($direct_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type, $info);   
    $set_unmapped_sth->execute($analysis_id, $ex_id, $acc, $reason_id{"NO_STABLE_ID"});
    push @xref_list, $xref_id;
  }
  $direct_unmapped_sth->finish;
  $set_unmapped_sth->finish;


  if(@xref_list){
    my $xref_dumped_sth = $self->xref->dbc->prepare("update xref set dumped = 4 where xref_id in (".join(", ",@xref_list).")");
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

  my $misc_unmapped_sth = $self->xref->dbc->prepare($sql);
  $misc_unmapped_sth->execute();
  $misc_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname);

  @xref_list = ();
  $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
  while($misc_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type, $info);   
    $set_unmapped_sth->execute($analysis_id, $ex_id, $acc, $reason_id{"NO_MAPPING"});
    push @xref_list, $xref_id;
  }
  $misc_unmapped_sth->finish;
  $set_unmapped_sth->finish;


  if(@xref_list){
    my $xref_dumped_sth = $self->xref->dbc->prepare("update xref set dumped = 5 where xref_id in (".join(", ",@xref_list).")");
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

  my $dep_unmapped_sth = $self->xref->dbc->prepare($sql);
  $dep_unmapped_sth->execute();
  my $parent;
  $dep_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname, \$parent);

  $set_unmapped_sth  =  $self->core->dbc->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id, parent ) values ('xref', ?, ?, ?, '".$reason_id{"MASTER_FAILED"}."', ?)");

  @xref_list = ();
  my $last_acc= 0;
  while($dep_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    if($last_acc ne $acc){
      $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label||$acc, $version, $desc, $type, $info);   
    }
    $last_acc = $acc;
    $set_unmapped_sth->execute($analysis_id, $ex_id, $acc, $parent);
    push @xref_list, $xref_id;
  }
  $dep_unmapped_sth->finish;
  $set_unmapped_sth->finish;

  if(@xref_list){
    my $xref_dumped_sth = $self->xref->dbc->prepare("update xref set dumped = 6 where xref_id in (".join(", ",@xref_list).")");
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
            ix.query_identity, ix.target_identity
      FROM source s, primary_xref px, xref x
        LEFT JOIN object_xref ox ON ox.xref_id = x.xref_id
        LEFT JOIN identity_xref ix ON ix.object_xref_id = ox.object_xref_id
      WHERE x.source_id = s.source_id
	  AND px.xref_id = x.xref_id
          AND x.dumped is null 
          AND x.info_type = 'SEQUENCE_MATCH'
          AND ox.ox_status != 'FAILED_PRIORITY'
          ORDER  BY x.xref_id
          
SEQ


  my $seq_unmapped_sth = $self->xref->dbc->prepare($sql);
  $seq_unmapped_sth->execute();
  my ($ensembl_object_type, $ensembl_id, $q_id, $t_id, $seq_type) ;
  $seq_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname, \$seq_type, \$ensembl_object_type, \$ensembl_id, \$q_id, \$t_id);

  my $set_unmapped_no_sth     = $self->core->dbc->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id, ensembl_object_type ) values ('xref', ?, ?, ?, '".$reason_id{"FAILED_MAP"}."', ?)");
  my $set_unmapped_failed_sth = $self->core->dbc->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id, query_score, target_score, ensembl_id, ensembl_object_type ) values ('xref', ?, ?, ?, '".$reason_id{"FAILED_THRESHOLD"}."',?,?,?,?)");


  @xref_list = ();
  my $last_xref = 0;
  while($seq_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    if($last_xref != $xref_id){
      $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type, $info);   
    }
    $last_xref = $xref_id;
    if(defined($ensembl_id)){
      $analysis_id= $analysis_ids{$ensembl_object_type};
      $set_unmapped_failed_sth->execute($analysis_id, $ex_id, $acc, $q_id, $t_id, $ensembl_id, $ensembl_object_type );
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
    my $xref_dumped_sth = $self->xref->dbc->prepare("update xref set dumped = 7 where xref_id in (".join(", ",@xref_list).")");
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

  
  my $wel_unmapped_sth = $self->xref->dbc->prepare($sql);
  $wel_unmapped_sth->execute();
  $wel_unmapped_sth->bind_columns(\$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname);

  $set_unmapped_sth  =  $self->core->dbc->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id) values ('xref', ?, ?, ?, '".$reason_id{"NO_MASTER"}."')");

  $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
  @xref_list = ();
  while($wel_unmapped_sth->fetch()){
    my $ex_id = $name_to_external_db_id{$dbname};
    $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type, $info);   
    $set_unmapped_sth->execute($analysis_id, $ex_id, $acc);
    push @xref_list, $xref_id;
  }
  $wel_unmapped_sth->finish;
  $set_unmapped_sth->finish;

  if(@xref_list){
    my $xref_dumped_sth = $self->xref->dbc->prepare("update xref set dumped = 8 where xref_id in (".join(", ",@xref_list).")");
    $xref_dumped_sth->execute(); 
    $xref_dumped_sth->finish;
  }


  $transaction_end_sth->execute();

  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('core_loaded',now())");
  $sth_stat->execute();
  $sth_stat->finish;



}


sub get_analysis{
  my $self = shift;
  

  my %typeToLogicName = ( 'Gene'        => 'XrefExonerateDNA',
                          'Transcript'  => 'XrefExonerateDNA',
                          'Translation' => 'XrefExonerateProtein' );

  my %analysis_id;

  foreach my $key (qw(Gene Transcript Translation)){
    
    my $logic_name = $typeToLogicName{$key};
    
    my $sth = $self->core->dbc->prepare("SELECT analysis_id FROM analysis WHERE logic_name='" . $logic_name ."'");
    
    $sth->execute();
    
    my $analysis_id;
    
    if (my @row = $sth->fetchrow_array()) {
      
      $analysis_id{$key} = $row[0];
      
    } else {
      
      print "No analysis with logic_name $logic_name found, creating ...\n" if ($self->verbose);
      $sth = $self->core->dbc->prepare("INSERT INTO analysis (logic_name, created) VALUES ('" . $logic_name. "',NOW())");
      # TODO - other fields in analysis table
      $sth->execute();
      $analysis_id{$key} = $sth->{'mysql_insertid'};
    }
    $sth->finish();
    
  }
  return %analysis_id;
  
}



1;

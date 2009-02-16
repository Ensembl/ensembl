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
  $sql = "select s.source_id, s.name from source s, xref x where x.source_id = s.source_id group by s.source_id"; # only get those of interest
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

  
  ######################################
  # For each external_db to be updated #
  # Delete the existing ones           # 
  ######################################


  $sth = $self->xref->dbc->prepare('select s.name, count(*) from xref x, object_xref ox, source s where ox.xref_id = x.xref_id  and x.source_id = s.source_id and ox_status = "DUMP_OUT"  group by s.name');
  $sth->execute();
  my $count;
  $sth->bind_columns(\$name,\$count);

  my $synonym_sth  =  $self->core->dbc->prepare('DELETE external_synonym FROM external_synonym, xref WHERE external_synonym.xref_id = xref.xref_id AND xref.external_db_id = ?');
  my $go_sth       =  $self->core->dbc->prepare('DELETE gx FROM xref x, object_xref ox LEFT JOIN go_xref gx ON ox.object_xref_id = gx.object_xref_id WHERE x.xref_id = ox.xref_id AND x.external_db_id = ? AND gx.linkage_type is not null');
  my $identity_sth =  $self->core->dbc->prepare('DELETE identity_xref FROM identity_xref, object_xref, xref WHERE identity_xref.object_xref_id = object_xref.object_xref_id AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = ?');
  my $object_sth   =  $self->core->dbc->prepare('DELETE object_xref FROM object_xref, xref WHERE object_xref.xref_id = xref.xref_id AND xref.external_db_id = ?');
#  my $dependent_sth = $self->core->dbc->prepare('DELETE dependent_xref FROM dependent_xref, xref  WHERE dependent_xref.dependent_xref_id = xref.xref_id and xref.external_db_id = ?');
  my $xref_sth     =  $self->core->dbc->prepare('DELETE FROM xref WHERE xref.external_db_id = ?');
  my $unmapped_sth =  $self->core->dbc->prepare('DELETE FROM unmapped_object WHERE type="xref" and external_db_id = ?');


  my $test =1;
  if(!$test){
  while($sth->fetch()){
    my $ex_id = $name_to_external_db_id{$name};

    print "Deleting data for $name from core before updating from new xref database\n";
    $synonym_sth->execute($ex_id);
    $go_sth->execute($ex_id);
    $identity_sth->execute($ex_id);
    $object_sth->execute($ex_id);  
#    $dependent_sth->execute($ex_id);
    $xref_sth->execute($ex_id);
    $unmapped_sth->execute($ex_id);
  }
  $sth->finish;
}
  $synonym_sth->finish;
  $go_sth->finish;  
  $identity_sth->finish;
  $object_sth->finish;  
#  $dependent_sth->finish;
  $xref_sth->finish;
  $unmapped_sth->finish; 

  ###############################################################
  ##### Create temp table dependent_xref (until schema changes) #
  ###############################################################

 
  $sql = (<<SQL);
  Create TABLE dependent_xref(
     object_xref_id         INT NOT NULL,
     master_xref_id         INT NOT NULL,
     dependent_xref_id      INT NOT NULL,

     PRIMARY KEY( master_xref_id ),
     KEY dependent ( dependent_xref_id )

   ) COLLATE=latin1_swedish_ci TYPE=MyISAM
SQL

  $sth = $self->core->dbc->prepare($sql);
  $sth->execute || die "Could not create temp table dependent_xref\n";
  $sth->finish;

  ##### Delete this ONLY after the gene/transcript display_xref and description calculations.


  ##########################################
  # Get the offsets for object_xref, xref  #
  ##########################################

  $sth = $self->core->dbc->prepare('select MAX(xref_id) from xref');
  my $xref_offset;
  $sth->execute;
  $sth->bind_columns(\$xref_offset);
  $sth->fetch();
  $sth->finish;

  $sth = $self->core->dbc->prepare('select MAX(object_xref_id) from object_xref');
  my $object_xref_offset;
  $sth->execute;
  $sth->bind_columns(\$object_xref_offset);
  $sth->fetch();
  $sth->finish;


  ####################
  # Get analysis id's 
  ####################

  my %analysis_id = $self->get_analysis(); # 


  print "xref offset is $xref_offset, object_xref offset is $object_xref_offset\n";

  #####################################
  # Now add the new ones              #
  #####################################

     ###########################
     # SQL to get data from xref
     ###########################

     my $direct_sth = $self->xref->dbc->prepare('select x.xref_id, x.accession, x.label, x.version, x.description, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type from xref x, object_xref ox  where ox.ox_status = "DUMP_OUT" and ox.xref_id = x.xref_id and x.source_id = ? and x.info_type = ? order by x.xref_id');
 
     my $dependent_sth = $self->xref->dbc->prepare('select  x.xref_id, x.accession, x.label, x.version, x.description, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, d.master_xref_id from xref x, object_xref ox,  dependent_xref d where ox.ox_status = "DUMP_OUT" and ox.xref_id = x.xref_id and d.dependent_xref_id = x.xref_id and x.source_id = ? and x.info_type = ? order by x.xref_id, ox.ensembl_id');

     $go_sth = $self->xref->dbc->prepare('select  x.xref_id, x.accession, x.label, x.version, x.description, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, d.master_xref_id, g.linkage_type from xref x, object_xref ox,  dependent_xref d, go_xref g where ox.ox_status = "DUMP_OUT" and  g.object_xref_id = ox.object_xref_id and x.xref_id = ox.xref_id and d.object_xref_id = ox.object_xref_id and x.source_id = ? and x.info_type = ? order by x.xref_id, ox.ensembl_id');

     my $seq_sth   =   $self->xref->dbc->prepare('select x.xref_id, x.accession, x.label, x.version, x.description, ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, i.query_identity, i.target_identity, i.hit_start, i.hit_end, i.translation_start, i.translation_end, i.cigar_line, i.score, i.evalue from xref x, object_xref ox, identity_xref i  where ox.ox_status = "DUMP_OUT" and i.object_xref_id = ox.object_xref_id and ox.xref_id = x.xref_id and x.source_id = ? and x.info_type = ? order by x.xref_id');

     ########################
     # SQL to add data to core
     #########################
 
     my $add_xref_sth           = $self->core->dbc->prepare('insert into xref (xref_id, external_db_id, dbprimary_acc, display_label, version, description, info_type) values (?, ?, ?, ?, ?, ?, ?)');
     my $add_object_xref_sth    = $self->core->dbc->prepare('insert into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id) values (?, ?, ?, ?)');
     my $add_identity_xref_sth  = $self->core->dbc->prepare('insert into identity_xref (object_xref_id, xref_identity, ensembl_identity, xref_start, xref_end, ensembl_start, ensembl_end, cigar_line, score, evalue, analysis_id) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
     my $add_go_xref_sth        = $self->core->dbc->prepare('insert into go_xref (object_xref_id, linkage_type) values (?, ?)');
     my $add_dependent_xref_sth = $self->core->dbc->prepare('insert into dependent_xref (object_xref_id, master_xref_id, dependent_xref_id) values (?, ?, ?)');
     my $add_syn_sth            = $self->core->dbc->prepare('insert into external_synonym (xref_id, synonym) values (?, ?)');

  $sth = $self->xref->dbc->prepare('select s.name, s.source_id, count(*), x.info_type from xref x, object_xref ox, source s where ox.xref_id = x.xref_id  and x.source_id = s.source_id and ox_status = "DUMP_OUT"  group by s.name, x.info_type');
  $sth->execute();
  my ($type, $source_id);
  $sth->bind_columns(\$name,\$source_id, \$count, \$type);
  while($sth->fetch()){
 
    my $ex_id = $name_to_external_db_id{$name};

    print "updating $name in core (for $type xrefs)\n";

    my @xref_list=();  # process at end. Add synonyms and set dumped = 1;

   
    # dump SEQUENCE_MATCH, DEPENDENT, DIRECT, COORDINATE_OVERLAP, INFERRED_PAIR, (MISC?? same as direct come from official naming)  

    ### If DIRECT ,         xref, object_xref,                  (order by xref_id)  # maybe linked to more than one?
    ### if INFERRED_PAIR    xref, object_xref
    ### if MISC             xref, object_xref 

    
    if($type eq "DIRECT" or $type eq "INFERRED_PAIR" or $type eq "MISC"){
      my $count = 0;
      $direct_sth->execute($source_id, $type);
      my ($xref_id, $acc, $label, $version, $desc, $object_xref_id, $ensembl_id, $ensembl_type); 
      $direct_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$object_xref_id, \$ensembl_id, \$ensembl_type);
      my $last_xref = 0;
      while($direct_sth->fetch){
        if($last_xref != $xref_id){
	  push @xref_list, $xref_id;
	  $count++;
	  $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type);
	  $last_xref = $xref_id;
        }
        $add_object_xref_sth->execute(($object_xref_id+$object_xref_offset), $ensembl_id, $ensembl_type, ($xref_id+$xref_offset));
      }  
      print "DIRECT $count\n";
    }
 
    ### If DEPENDENT,       xref, object_xref , dependent_xref  (order by xref_id)  # maybe linked to more than one?
 
   elsif($type eq "DEPENDENT"){
     if($name eq "GO"){
       my $count = 0;
       $go_sth->execute($source_id, $type);
       my ($xref_id, $acc, $label, $version, $desc, $object_xref_id, $ensembl_id, $ensembl_type, $master_xref_id, $linkage_type); 
       $go_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$object_xref_id, \$ensembl_id, \$ensembl_type, \$master_xref_id, \$linkage_type);
       my $last_xref = 0;
       while($go_sth->fetch){
	 if($last_xref != $xref_id){
	   push @xref_list, $xref_id;
	   $count++;
	   $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type);
	   $last_xref = $xref_id;
	 }
	 $add_dependent_xref_sth->execute(($object_xref_id+$object_xref_offset), ($xref_id+$xref_offset), ($master_xref_id+$xref_offset) );
	 $add_object_xref_sth->execute( ($object_xref_id+$object_xref_offset), $ensembl_id, $ensembl_type, ($xref_id+$xref_offset) );
	 $add_go_xref_sth->execute( ($object_xref_offset+$object_xref_id), $linkage_type);
       }       
       print "GO $count\n";     
     }
    else{
      my $count = 0;
      $dependent_sth->execute($source_id, $type);
      my ($xref_id, $acc, $label, $version, $desc, $object_xref_id, $ensembl_id, $ensembl_type, $master_xref_id); 
      $dependent_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$object_xref_id, \$ensembl_id, \$ensembl_type, \$master_xref_id);
      my $last_xref = 0;
      my $last_ensembl = 0;
      while($dependent_sth->fetch){
        if($last_xref != $xref_id){
	  push @xref_list, $xref_id;
	  $count++;
	  $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label || $acc, $version, $desc, $type);
	  $last_xref = $xref_id;
        }
	if($last_xref != $xref_id or $last_ensembl != $ensembl_id){
	  $add_object_xref_sth->execute(($object_xref_id+$object_xref_offset), $ensembl_id, $ensembl_type, ($xref_id+$xref_offset));
	  $add_dependent_xref_sth->execute(($object_xref_id+$object_xref_offset), ($xref_id+$xref_offset), ($master_xref_id+$xref_offset) );	}
	$last_ensembl = $ensembl_id;
      }  
      print "DEP $count\n";
    }
   }
   ### If SEQUENCE_MATCH   xref, object_xref,  identity_xref   (order by xref_id)  # maybe linked to more than one?

    elsif($type eq "SEQUENCE_MATCH"){
      my $count = 0;
      $seq_sth->execute($source_id, $type);
      my ($xref_id, $acc, $label, $version, $desc, $object_xref_id, $ensembl_id, $ensembl_type); 
      my ( $query_identity, $target_identity, $hit_start, $hit_end, $translation_start, $translation_end, $cigar_line, $score, $evalue);
      $seq_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$object_xref_id, \$ensembl_id, \$ensembl_type,
			     \$query_identity, \$target_identity, \$hit_start, \$hit_end, \$translation_start, \$translation_end, \$cigar_line, \$score, \$evalue);
      my $last_xref = 0;
      while($seq_sth->fetch){
        if($last_xref != $xref_id){
	  push @xref_list, $xref_id;
	  $count++;
	  $add_xref_sth->execute(($xref_id+$xref_offset), $ex_id, $acc, $label, $version, $desc, $type);
	  $last_xref = $xref_id;
        }
        $add_object_xref_sth->execute(($object_xref_id+$object_xref_offset), $ensembl_id, $ensembl_type, ($xref_id+$xref_offset));
	$add_identity_xref_sth->execute( ($object_xref_id+$object_xref_offset), $query_identity, $target_identity, $hit_start, $hit_end, 
					 $translation_start, $translation_end, $cigar_line, $score, $evalue, $analysis_id{$ensembl_type});  
      }  
      print "SEQ $count\n";
    }
    else{
      print "ARSE what type is $type\n";
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

      my $xref_dumped_sth = $self->xref->dbc->prepare("update xref set dumped = 1 where xref_id in (".join(", ",@xref_list).")");
      $xref_dumped_sth->execute(); 
      $xref_dumped_sth->finish;
    }	



    # if its a priority xref :-
       # write unmapped xrefs    


    # else not priority xref
       # write unmapped xrefs
 


  }
  $sth->finish;


  # remove after testing
  $sth = $self->core->dbc->prepare("drop table dependent_xref");
  $sth->execute || die "Could not drop temp table dependent_xref\n";
  $sth->finish;  

}


sub get_analysis{
  my $self = shift;
  

  my %typeToLogicName = ( 'Transcript' => 'XrefExonerateDNA',
                          'Translation' => 'XrefExonerateProtein' );

  my %analysis_id;

  foreach my $key (qw(Transcript Translation)){
    
    my $logic_name = $typeToLogicName{$key};
    
    my $sth = $self->core->dbc->prepare("SELECT analysis_id FROM analysis WHERE logic_name='" . $logic_name ."'");
    
    $sth->execute();
    
    my $analysis_id;
    
    if (my @row = $sth->fetchrow_array()) {
      
      $analysis_id{$key} = $row[0];
      
    } else {
      
      print "No analysis with logic_name $logic_name found, creating ...\n";
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

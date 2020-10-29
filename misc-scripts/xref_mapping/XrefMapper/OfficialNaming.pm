=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

package XrefMapper::OfficialNaming;

use strict;
use warnings;
use Carp;
use Cwd;
use DBI;
use File::Basename;

use vars '@ISA';
use base qw( XrefMapper::BasicMapper);
#@ISA = qw{ XrefMapper::BasicMapper };


###############################################################################################
# Run the offical naming code.
#
#   At present this is done for the following species ONLY :-
#      ZebraFish (ZFIN_ID),
#      Human (HGNC)
#      Mouse (MGI)
#      Rat (RGD)
#      Pig (PIGGY)
#         There is currently no official domain source for pig, but it has manual annotation
#         We use PIGGY as a fake official naming source
#
#  1) So we find the best official name for each gene
#     order for this is:-
#               i)   official domain name source (HGNC, MGI, ZFIN_ID, RGD)
#               ii)  RFAM
#               iii) miRBase
#               iv)   EntrezGene
#               v) Clone name
#
#      NOTE: for "i)" above, if more than one exists we find the "best" one if possible
#            and remove the other ones. If there is more than one "best" we keep all and
#            just choose the first one for the name
#
#            To find the "best" one we use the priority. 
#            Priority should be set correctly in the xref_config.ini file to use
#            first any names coming from the official naming source
#
#      Set this as the display_xref for the gene.
#
#  2) Foreach Transcript of that gene 
#
#     we assign a transcript extension (splice number?)
#     This is just a counter starting at 201 which is incremented each time
#     We add this to the name to get a "XXX_trans_name"xref  where XXX is the 
#     type of source used to get the name. This is then added as an xref and 
#     is set to the display_xref for that transcript.
#
##############################################################################################


####################################
# Create OfficialNaming object
# Get some info from the BasicMapper
####################################
sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class; 
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->get_official_name($mapper->get_official_name);
  return $self;
}


##################################################
# This will be the offical database name
# HGNC, MGI, ZFIN_ID or PIGGY, comes from BasicMapper
#################################################
sub get_official_name {
 my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_official_name} = $arg );
  return $self->{_official_name};
}



##################################################
# This is the main subroutine that does everything
##################################################
sub run {
  my $self = shift;
  my $species_id = shift;

  my $dbname = $self->get_official_name();
  my $dbi = $self->xref->dbc;

  ###########################################################
  # If no offical name then we do not want to go any further
  # Just set status to official_naming_done and return
  ###########################################################
  if(!defined($dbname)){
    $self->update_process_status("official_naming_done");
    return;
  }
  $species_id = $self->get_id_from_species_name($self->core->species) unless defined $species_id;
  $self->species_id($species_id);


  ###########################################################
  # If there are any official names on transcripts or translations
  # move them onto gene level
  #
  # This is done for 2 reasons
  #  1) to make the code the same as HGNC is on a gene
  #     and it makes it easier to find.
  #  2) Later on these are copied to the canonical transcripts
  #     from the genes so move them now.
  ###########################################################

  if($dbname eq "MGI"){ # Copy MGI to Genes
    $self->biomart_fix("MGI","Translation","Gene");
    $self->biomart_fix("MGI","Transcript","Gene");
  }
  if($dbname eq "ZFIN_ID"){ # Copy ZFIN_ID to Genes
    $self->biomart_fix("ZFIN_ID","Translation","Gene");
    $self->biomart_fix("ZFIN_ID","Transcript","Gene");
  }
  if($dbname eq "RGD"){ # Copy RGD to Genes
    $self->biomart_fix("RGD","Translation","Gene");
    $self->biomart_fix("RGD","Transcript","Gene");
  }



  ######################################################
  # Get the current max values for xref and object_xref
  ######################################################
  my ($max_object_xref_id, $max_xref_id) = $self->find_max_ids($dbi);

  my %display_label_to_desc;
  $self->get_display_label_data(\%display_label_to_desc, $dbi);

  my %synonym;
  $self->get_synonyms(\%synonym, $dbi);


  # get the officail naming external_sources
  my $dbname_to_source_id = $self->get_new_dbname_sources($dbi); # reference to hash

  ###########################
  # Delete the old ones.
  ###########################
  $self->delete_old_data($dbname_to_source_id, $dbi);

  $self->reset_display_xrefs($dbi);

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbconn => $self->core->dbc);
  my $ga = $db->get_GeneAdaptor();

  my %gene_to_transcripts;
  my %gene_id_to_stable_id;
  my %tran_id_to_stable_id;

  my $sql =(<<'SQ0');
SELECT gtt.gene_id, gtt.transcript_id, gsi.stable_id, tsi.stable_id 
  FROM gene_transcript_translation gtt, gene_stable_id gsi, transcript_stable_id tsi 
    WHERE gtt.gene_id = gsi.internal_id AND
          gtt.transcript_id = tsi.internal_id
    ORDER BY gsi.stable_id, tsi.stable_id
SQ0

  my $sth = $dbi->prepare($sql);

  $sth->execute;
  my ($gene_id, $tran_id, $gsi, $tsi);
  $sth->bind_columns(\$gene_id, \$tran_id, \$gsi, \$tsi);
  my @sorted_gene_ids;
  while ($sth->fetch){
    if(!defined($gene_to_transcripts{$gene_id})){
      push @sorted_gene_ids, $gene_id;
    }
    push @{$gene_to_transcripts{$gene_id}}, $tran_id;
    $gene_id_to_stable_id{$gene_id} = $gsi; 
    $tran_id_to_stable_id{$tran_id} = $tsi; 
  }

  my $dbentrie_sth = $self->get_dbentrie_sth($dbi);
  my $ins_xref_sth = $self->get_ins_xref_sth($dbi);
  my $ins_dep_ix_sth = $self->get_ins_dep_ix_sth($dbi);
  my $ins_object_xref_sth =  $self->get_ins_object_xref_sth($dbi);
  my $set_gene_display_xref_sth = $self->get_set_gene_display_xref_sth($dbi);

  my %xref_added; # store those added  $xref_added{$accession:$source_id} = $xref_id;
  my %seen_gene;

  my %official_name_used;

  my $ignore_sql =<<IEG;
  SELECT DISTINCT ox.object_xref_id
    FROM object_xref ox, dependent_xref dx,
       xref xmas, xref xdep,
       source smas, source sdep
    WHERE ox.xref_id = dx.dependent_xref_id AND
          dx.dependent_xref_id = xdep.xref_id AND
          dx.master_xref_id = xmas.xref_id AND
          xmas.source_id = smas.source_id AND
          xdep.source_id = sdep.source_id AND
          smas.name like "Refseq%predicted" AND
          sdep.name like "EntrezGene" AND
          ox.ox_status = "DUMP_OUT"
IEG

  my %ignore_object;
  my $ignore_sth = $dbi->prepare($ignore_sql);
  $ignore_sth->execute();
  my ($ignore_object_xref_id);
  $ignore_sth->bind_columns(\$ignore_object_xref_id);
  while($ignore_sth->fetch()){
    $ignore_object{$ignore_object_xref_id} = 1;
  }
  $ignore_sth->finish;

  while ( my $gene_id = shift @sorted_gene_ids){

    my $tran_source = $dbname;

    # symbols to set when found.
    my $gene_symbol = undef;
    my $gene_symbol_xref_id = undef;
    my $is_lrg = 0;

    ################################
    # Get offical name if it has one
    ################################
    ($gene_symbol, $gene_symbol_xref_id) = 
      $self->get_official_domain_name({gene_id       => $gene_id, 
				       gene_to_tran  => \%gene_to_transcripts,
				       gene_id_to_stable_id => \%gene_id_to_stable_id,
                                       official_name_used => \%official_name_used,
                                       dbi            => $dbi
				      });

    if (defined($gene_symbol_xref_id)) {
	$official_name_used{$gene_symbol_xref_id} = 1;
    }

    ############################################
    # If not found see if there is an LRG entry
    ############################################
    if(!defined($gene_symbol)){ # look for LRG
      ($gene_symbol, $gene_symbol_xref_id, $is_lrg) = $self->find_lrg_hgnc($gene_id, $dbi);
    }

    ####################################################
    # If not found look for other valid database sources
    # These are RFAM and miRBase, as well as EntrezGene
    ####################################################
    if(!defined($gene_symbol)){ 
      ($gene_symbol, $gene_symbol_xref_id) = 
	  $self->find_from_other_sources(\%ignore_object, 
                                     {gene_id       => $gene_id, 
					label_to_desc => \%display_label_to_desc,
                                      dbi           => $dbi,
					tran_source   => \$tran_source});
    }

    if(defined($gene_symbol)){
      my $desc = $display_label_to_desc{$gene_symbol};
      $set_gene_display_xref_sth->execute($gene_symbol_xref_id, $gene_id);

      if (!$is_lrg) {
        $self->set_transcript_display_xrefs({ max_xref         => \$max_xref_id, 
                                            max_object       => \$max_object_xref_id,
                                            gene_id          =>  $gene_id,
                                            gene_id_to_stable_id => \%gene_id_to_stable_id,
					    gene_symbol      => $gene_symbol,
					    desc             => $desc, 
                                            dbi              => $dbi,
                                            source_id        => $dbname_to_source_id->{$tran_source."_trans_name"}, 
                                            xref_added       => \%xref_added, 
                                            seen_gene        => \%seen_gene, 
                                            gene_to_tran     => \%gene_to_transcripts, 
                                            tran_source      => $tran_source,
					   });
      }
    }

  } # for each gene

  $self->update_process_status('official_naming_done');
  return;
}




####################################################################
# Get offical name if it has one
#
# Search gene for dbname entries.
# dbname (HGNC||MGI||ZFIN_ID|RGD) dependent on species
#
# Find the "best" one
# Remove the lesser ones (set status to MULTI_DELETE for object_xref)
#
# return the gene_symbol and xref_id of the best one
######################################################################

sub get_official_domain_name{
  my ($self, $arg_ref) = @_;

  my $gene_id              = $arg_ref->{gene_id};
  my $gene_id_to_stable_id = $arg_ref->{gene_id_to_stable_id};
  my $gene_to_transcripts  = $arg_ref->{gene_to_tran};
  my $official_name_used   = $arg_ref->{official_name_used};
  my $dbi                  = $arg_ref->{dbi};


  my $dbname = $self->get_official_name();
  my $gene_symbol = undef;
  my $gene_symbol_xref_id = undef;


  my $dbentrie_sth = $self->get_dbentrie_sth($dbi);

  my %ODN=();
  my %xref_id_to_display;

  $dbentrie_sth->execute($dbname, $gene_id, "Gene");
  my ($display, $xref_id, $object_xref_id, $level); 
  $dbentrie_sth->bind_columns(\$display, \$xref_id, \$object_xref_id, \$level);
  my $best_level=999;
  
  my $count = 0;
  my @list=();
  my @list_ox=();

  while($dbentrie_sth->fetch){

    push @list, $xref_id;
    push @list_ox, $object_xref_id;
    $count++;
    $xref_id_to_display{$xref_id} = $display;
    if($level < $best_level){
      %ODN = ();
      $ODN{$xref_id} = 1;
      $best_level = $level;
    }
    elsif($level == $best_level){
      $ODN{$xref_id} = 1;
    }
  }

  if(($count > 1) and (scalar(keys %ODN) == 1)){ # found one that is "best" so set it and remove others
    print "For gene ".$gene_id_to_stable_id->{$gene_id}." we have mutiple ".$dbname."'s\n";
    ($gene_symbol, $gene_symbol_xref_id) = $self->set_the_best_odns(\%ODN, \@list, \@list_ox, \%xref_id_to_display, $dbi);
    if(defined($gene_symbol)){
      return $gene_symbol, $gene_symbol_xref_id;
    }
  }

  if(scalar(keys %ODN) == 1){  # one hgnc to this gene - perfect case :-)
      return $xref_id_to_display{(keys %ODN)[0]}, (keys %ODN)[0];
  }
  if(scalar(keys %ODN) > 1){ 
    
    #if we have  more than 1 xref, fail xrefs with worse % identity if we can (query or target identity whichever is greater)
    my $identity_sth = $self->get_best_identity_sth($dbi);     
    $identity_sth->execute($dbname, $gene_id, "Gene");
    my ($xref_id, $best_identity); 
    $identity_sth->bind_columns(\$xref_id, \$best_identity);
    my $temp_best_identity = 0;
    my %best_ids = ();

    while($identity_sth->fetch){
	 
	if($best_identity > $temp_best_identity){
	     %best_ids = ();
	     $best_ids{$xref_id} = 1;
	     $temp_best_identity = $best_identity;
	}
	elsif($best_identity == $temp_best_identity){
	     $best_ids{$xref_id} = 1;
	} 
	else {
	     last;
	}
    }

    my %best_list;
    foreach my $xref_id (keys %ODN){
      $best_list{$xref_id_to_display{$xref_id}} = 1;
    }

    # check if we were able to reduce the number of xrefs based on % identity 
    if ( scalar(keys %best_ids) > 0 && scalar(keys %best_ids) < scalar(keys %ODN) ) {
	  %ODN = %best_ids;
	  print "For gene ".$gene_id_to_stable_id->{$gene_id}." we have mutiple ".$dbname."'s\n";
	  #set statuses for xrefs with worse % identity to MULTI_DELETE
	  ($gene_symbol, $gene_symbol_xref_id) = $self->set_the_best_odns(\%ODN, \@list, \@list_ox, \%xref_id_to_display, $dbi);
	  if( defined($gene_symbol) && scalar(keys %ODN == 1) ){
	      return $gene_symbol, $gene_symbol_xref_id;
	  } 
    }
    
    # take the name which hasn't been already assigned to another gene, if possible
    
    my $xref_not_used;
    foreach my $x (keys %ODN){
	  if (!defined($official_name_used->{$x}) ) {
	      $xref_not_used = $x;
	  }
    }
    if ($xref_not_used) {
	  foreach my $x (keys %ODN){
	      print "\t".$xref_id_to_display{$x};
	      if ($x == $xref_not_used) {
		  print "    chosen\n";
		  $gene_symbol =  $xref_id_to_display{$x};
		  $gene_symbol_xref_id =  $x;		  
	      } else {
		  print "  (left as $dbname reference but not gene symbol)\n";
	      }
	  }

    } else {

	  my $i=0;
	  foreach my $x (keys %ODN){
	      print "\t".$xref_id_to_display{$x};
	      if(!$i){
		  print "  (chosen as first)\n";
		  $gene_symbol =  $xref_id_to_display{$x};
		  $gene_symbol_xref_id =  $x;
	      }
	      else{
		  print "  (left as $dbname reference but not gene symbol)\n";
	      }
	      $i++;
	  }


    }
  }
  return ($gene_symbol, $gene_symbol_xref_id);
}


###########################################################
# Set the transcript display xrefs
#
# Use the gene symbol to create a transcript display xref
# Add 201 and increment.
###########################################################
sub set_transcript_display_xrefs{
  my ($self, $arg_ref) = @_;

  my $max_xref_id =         $arg_ref->{max_xref};
  my $max_object_xref_id =  $arg_ref->{max_object};
  my $gene_id =             $arg_ref->{gene_id};
  my $gene_symbol =         $arg_ref->{gene_symbol};
  my $desc =                $arg_ref->{desc};
  my $source_id =           $arg_ref->{source_id};
  my $xref_added =          $arg_ref->{xref_added};
  my $seen_gene =           $arg_ref->{seen_gene};
  my $gene_to_transcripts = $arg_ref->{gene_to_tran};
  my $tran_source         = $arg_ref->{tran_source};
  my $gene_id_to_stable_id = $arg_ref->{gene_id_to_stable_id};
  my $dbi                  = $arg_ref->{dbi};


  # statement handles needed
  my $ins_xref_sth =              $self->get_ins_xref_sth($dbi);
  my $ins_dep_ix_sth =            $self->get_ins_dep_ix_sth($dbi);
  my $set_tran_display_xref_sth = $self->get_set_transcript_display_xref_sth($dbi);
  my $ins_object_xref_sth =       $self->get_ins_object_xref_sth($dbi);

  if ($gene_id_to_stable_id->{$gene_id} =~ /LRG/) { return; }

  my $ext = 201;
  if(defined($seen_gene->{$gene_symbol})){
    $ext = $seen_gene->{$gene_symbol};
  }

  foreach my $tran_id ( @{$gene_to_transcripts->{$gene_id}} ){
    my $id = $gene_symbol."-".$ext;
    if(!defined($source_id)){
      croak "id = $id\n but NO source_id for this entry for $tran_source???\n";
    }
    if(!defined($xref_added->{$id.":".$source_id})){
      $$max_xref_id++;
      $ins_xref_sth->execute($$max_xref_id, $source_id, $id, $id, "", $desc);
      $xref_added->{$id.":".$source_id} = $$max_xref_id;
    }
    $set_tran_display_xref_sth->execute($xref_added->{$id.":".$source_id}, $tran_id);
    $$max_object_xref_id++;
    $ins_object_xref_sth->execute($$max_object_xref_id, $tran_id, 'Transcript', $xref_added->{$id.":".$source_id},undef);
    $ins_dep_ix_sth->execute($$max_object_xref_id, 100, 100);
    $ext++;
  }
  $seen_gene->{$gene_symbol} = $ext;
  return;
}


#################################################
# Get statement handle to retrieve what xrefs
# are attached to a specific ensembl_id and type
# for a particular source name
#################################################
sub get_dbentrie_sth{
  my $self = shift;
  my $dbi = shift;


  my $sql =(<<"SQ1");
SELECT x.label, x.xref_id, ox.object_xref_id, s.priority 
  FROM xref x, object_xref ox, source s
    WHERE x.xref_id = ox.xref_id AND
          x.source_id = s.source_id AND
          s.name = ? AND
          ox.ox_status = 'DUMP_OUT' AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = ?
SQ1
  my $sth = $dbi->prepare($sql);
  return $sth;
}

#################################################
# Get statement handle to retrieve what xrefs
# are attached to a specific ensembl_id and type
# for a particular source name with description
#################################################
sub get_dbentrie_with_desc_sth{
  my $self = shift;
  my $dbi = shift;


  my $sql =(<<"SQD");
SELECT x.label, x.xref_id, ox.object_xref_id, s.priority, x.description 
  FROM xref x, object_xref ox, source s
    WHERE x.xref_id = ox.xref_id AND
          x.source_id = s.source_id AND
          s.name = ? AND
          ox.ox_status = 'DUMP_OUT' AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = ?
SQD
  my $sth = $dbi->prepare($sql);
  return $sth;
}

#################################################
# Get statement handle to retrieve average of query
# and target identity for xrefs
#################################################
sub get_best_identity_sth{
  my $self = shift;
  my $dbi = shift;

  my $sql =(<<"SQD");
SELECT x.xref_id, CASE WHEN ix.query_identity >= ix.target_identity 
THEN ix.query_identity ELSE ix.target_identity END as best_identity 
FROM xref x, object_xref ox, identity_xref ix, source s 
WHERE x.xref_id = ox.xref_id AND x.source_id = s.source_id 
 AND ox.object_xref_id = ix.object_xref_id AND s.name = ? 
 AND ox.ox_status = 'DUMP_OUT' AND ox.ensembl_id = ? 
 AND ox.ensembl_object_type = ? order by best_identity DESC
SQD
  my $sth = $dbi->prepare($sql);
  return $sth;
}


#################################################
# Get statement handle to set the display xref
# for a transcript in the xref database.
# Stored in the transcript_stable_id table.
#################################################
sub get_set_transcript_display_xref_sth {
  my $self = shift;
  my $dbi = shift;
  my $sth = $dbi->prepare('UPDATE transcript_stable_id SET display_xref_id =? where internal_id = ?');
  return $sth;
}


#################################################
# Get statement handle to set the display xref
# for a gene in the xref database.
# Stored in the gene_stable_id table.
#################################################
sub get_set_gene_display_xref_sth {
  my $self = shift;
  my $dbi = shift;
  my $sth = $dbi->prepare('UPDATE gene_stable_id SET display_xref_id =? where internal_id = ?');
  return $sth;
}


###############################################
# Get statement handle to insert an xref
############################################### 
sub get_ins_xref_sth{
  my $self= shift;
  my $dbi = shift;

  my $sql = "insert ignore into xref (xref_id, source_id, accession, label, version, species_id, info_type, info_text, description) values (?, ?, ?, ?,  0, ".$self->species_id.", 'MISC', ?, ? )";
  my $sth = $dbi->prepare($sql); 
  return $sth;
}


#################################################
# Get statement handle to insert an identity xref
#################################################
sub get_ins_dep_ix_sth{
  my $self= shift;
  my $dbi = shift;

  my $sql = "insert into identity_xref (object_xref_id, query_identity, target_identity) values(?, ?, ?)";
  my $sth = $dbi->prepare($sql); 
  return $sth;
}

###############################################
# Get statement handle to insert an object_xref
############################################### 
sub get_ins_object_xref_sth{
  my $self= shift;
  my $dbi = shift;

  my $sql = "insert into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id, linkage_type, ox_status, unused_priority) values (?, ?, ?, ?, 'MISC', 'DUMP_OUT', ?)";
  my $sth = $dbi->prepare($sql); 
  return $sth;
}



sub find_max_ids{
  my $self = shift;
  my $dbi = shift;

  my ($max_object_xref_id, $max_object_xref_id2, $max_xref_id);

  my $sth = $dbi->prepare("SELECT MAX(object_xref_id) FROM object_xref");
  $sth->execute();
  $sth->bind_columns(\$max_object_xref_id);
  $sth->fetch;

  $sth = $dbi->prepare("SELECT MAX(object_xref_id) FROM identity_xref");
  $sth->execute();
  $sth->bind_columns(\$max_object_xref_id2);
  $sth->fetch;

  
  
  $sth = $dbi->prepare("SELECT MAX(xref_id) FROM xref");
  $sth->execute();
  $sth->bind_columns(\$max_xref_id);
  $sth->fetch;

  print "MAX xref_id = $max_xref_id MAX object_xref_id = $max_object_xref_id, max_object_xref from identity_xref = $max_object_xref_id2\n";
  return $max_object_xref_id, $max_xref_id;
}

sub get_synonyms{
  my ($self, $synonym, $dbi) = @_;

  my $dbname = $self->get_official_name();

  my $syn_sql = (<<"SYN");
SELECT es.synonym, x.label 
  FROM synonym es, xref x, source s 
    WHERE x.xref_id = es.xref_id AND
          x.source_id = s.source_id AND
           s.name = '$dbname'
SYN

  my $sth = $dbi->prepare($syn_sql);
  $sth->execute();
  my ($syn, $name);
  $sth->bind_columns(\$syn,\$name);
  while($sth->fetch){
    $synonym->{$syn} = $name;
  }
  $sth->finish;
  return;
}

sub get_display_label_data{
#  my ($self, $label_to_id, $label_to_desc) = @_;
  my ($self, $label_to_desc, $dbi) = @_;

  my $dbname = $self->get_official_name();

  my $gd1_sql = (<<"GD1");
SELECT x.accession, sy.synonym, x.description
  FROM synonym sy, xref x, source so
    WHERE x.xref_id = sy.xref_id AND
          so.source_id = x.source_id AND
           so.name like '$dbname'
GD1

  my $gd1_sth = $dbi->prepare($gd1_sql);

  $gd1_sth->execute();
  my ($display_label, $acc, $syn, $desc);
  $gd1_sth->bind_columns(\$acc,\$display_label, \$desc);
  while($gd1_sth->fetch){
#    $label_to_id->{$display_label} = $acc;
    $label_to_desc->{$display_label} = $desc;
  }
  $gd1_sth->finish;



  # get label to id from xref database to start with.
  my $gd2_sql = (<<"GD2");
SELECT x.accession, x.label, x.description
   FROM xref x, source s
     WHERE s.source_id = x.source_id AND
            s.name like '$dbname'
GD2
  
  my $gd2_sth = $dbi->prepare($gd2_sql);
  
  $gd2_sth->execute();
  $gd2_sth->bind_columns(\$acc,\$display_label, \$desc);
  while($gd2_sth->fetch){
#    $label_to_id->{$display_label} = $acc;
    if(!defined($desc)){
      warn "undef desc for $display_label\n";
    }
    else{
      $label_to_desc->{$display_label} = $desc;
    }
  }
  $gd2_sth->finish;
  return;
}

sub get_other_name_hash{
  my $self = shift;

  if(!defined($self->{'_other_name'})){
    my %hash;
    $self->{'_other_name'} = \%hash;
  }
  return  $self->{'_other_name'};
}




sub find_from_other_sources{
  my ($self, $ignore_object, $ref_args) = @_;
  my $tran_source           = $ref_args->{tran_source};
  my $gene_id               = $ref_args->{gene_id};
  my $display_label_to_desc = $ref_args->{label_to_desc}; 
  my $dbi                   = $ref_args->{dbi};
  my %ignore_object = %{$ignore_object};

  my ($gene_symbol, $gene_symbol_xref_id);
  my $dbentrie_sth = $self->get_dbentrie_with_desc_sth($dbi);
  my $other_name_num = $self->get_other_name_hash();

  my ($display, $xref_id, $object_xref_id, $level, $desc);
  my %found_gene;
  foreach my $ext_db_name (qw(miRBase RFAM EntrezGene)){
    $dbentrie_sth->execute($ext_db_name, $gene_id, "Gene");
    $dbentrie_sth->bind_columns(\$display, \$xref_id, \$object_xref_id, \$level, \$desc);
    while($dbentrie_sth->fetch){
      if (defined $found_gene{$gene_id}) {
        last;
      }
      if ($display =~ /^LOC/ || $display =~ /^SSC/) {
        next;
      }
      if (defined $ignore_object{$object_xref_id}) {
        next;
      }
      $gene_symbol = $display;
      $gene_symbol_xref_id = $xref_id;
      $$tran_source = $ext_db_name;
      $display_label_to_desc->{$display} = $desc;
      if(defined($other_name_num->{$gene_symbol})){
	$other_name_num->{$gene_symbol}++;
      }
      else{
	$other_name_num->{$gene_symbol} = 1;
      }
      if ($ext_db_name eq 'miRBase' || $ext_db_name eq 'RFAM') {
        $gene_symbol .= ".".$other_name_num->{$gene_symbol};
      }
      $found_gene{$gene_id} = 1;
      next;
    }
  }  
  return ($gene_symbol, $gene_symbol_xref_id);
}


#
# We do not delete this but set the status to "MULTI_DELETE"
#

sub get_delete_odn_sth{
  my ($self, $dbi) = @_;

  my $sth = $dbi->prepare('UPDATE object_xref SET ox_status = "MULTI_DELETE" where object_xref_id = ?');
  return $sth;
}

sub set_the_best_odns{
  my ($self, $odn, $ref_list, $ref_list_ox, $ref_xref_id_to_display, $dbi) = @_;

  my $delete_odn_sth = $self->get_delete_odn_sth($dbi);
  my %ODN = %$odn;

  my $gene_symbol = undef;
  my $gene_symbol_xref_id = undef;
  my $i=0;
  while ($i < scalar(@{$ref_list})){
    my $x = $ref_list->[$i];
    if(!exists($ODN{$x})){
      print "\tremoving ".$ref_xref_id_to_display->{$x}." from gene\n";
      #remove object xref....
      $delete_odn_sth->execute($ref_list_ox->[$i])|| 
	croak "Could not set staus to MULTI_DELETE for object_xref ".$ref_list_ox->[$i]."\n";
    }
    else{
      print "\tKeeping the best one ".$ref_xref_id_to_display->{$x}."\n";
      $gene_symbol = $ref_xref_id_to_display->{$x};
      $gene_symbol_xref_id = $x;
    }
    $i++;
  }
  return ($gene_symbol, $gene_symbol_xref_id);
}

########################## START LRG BIT ######################################################

sub get_lrg_find_sth{
  my $self = shift;
  my $dbi = shift;
  
  my $sql=(<<'SQ2');
SELECT x.label, x.xref_id, ox.object_xref_id, s.priority 
  FROM xref x, object_xref ox, source s 
    WHERE x.xref_id = ox.xref_id AND
          x.source_id = s.source_id AND 
          s.name = ? AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = ?
SQ2
  my $sth = $dbi->prepare($sql);
  return $sth;
}


sub get_lrg_set_status_sth{
  my $self = shift;
  my $dbi = shift;
  
  my $sth = $dbi->prepare("update object_xref set ox_status = 'NO_DISPLAY' where object_xref_id = ?");
  return $sth;
}

sub get_lrg_to_hgnc_sth{
  my $self = shift;
  my $dbi = shift;
  
  my $sql=(<<'SQ4');
SELECT x.xref_id, s.priority 
  FROM xref x,source s, object_xref ox
    WHERE x.xref_id = ox.xref_id AND
          x.source_id = s.source_id AND
          x.label = ? AND
          s.name = ? AND
          ox.ox_status = 'DUMP_OUT'
    ORDER BY s.priority
SQ4
  my $sth = $dbi->prepare($sql);
  return $sth;
}


sub find_lrg_hgnc{
  my ($self, $gene_id, $dbi) =@_;
  my $gene_symbol;
  my $gene_symbol_xref_id;
  my $is_lrg = 0;

  my $lrg_find_sth = $self->get_lrg_find_sth($dbi);
  my $lrg_set_status_sth = $self->get_lrg_set_status_sth($dbi);
  my $lrg_to_hgnc_sth = $self->get_lrg_to_hgnc_sth($dbi);

  # look for LRG_HGNC_notransfer, if found then find HGNC equiv and set to this
  #      print "LRG FOUND with no HGNC, should have gotten this via the alt allele table?? gene_id = $gene_id\n";
  $lrg_find_sth->execute("LRG_HGNC_notransfer", $gene_id, "Gene");
  my ($display, $xref_id, $object_xref_id, $level);
  $lrg_find_sth->bind_columns(\$display, \$xref_id, \$object_xref_id, \$level);
  while($lrg_find_sth->fetch){
    $lrg_set_status_sth->execute($object_xref_id); # set oc_status to no _display as we do not want this transferred, 
                                                   # just the equivalent hgnc
    my $new_xref_id  = undef;
    my $pp;
    $lrg_to_hgnc_sth->execute($display,"HGNC");
	$lrg_to_hgnc_sth->bind_columns(\$new_xref_id,\$pp);
    $lrg_to_hgnc_sth->fetch;
    if(defined($new_xref_id)){
      $gene_symbol = $display;
      $gene_symbol_xref_id = $new_xref_id;
      $is_lrg = 1;
    }
  }
  return ($gene_symbol, $gene_symbol_xref_id, $is_lrg);
}

#############################END LRG BIT ################################################

#
# These are the ones added by official naming and hence
# Need to be removed incase they still exist from a previous run
#
sub get_new_dbname_sources{
  my $self = shift;
  my $dbi = shift;

  my %dbname_to_source_id;

  my $dbname = $self->get_official_name();

  my @list = qw(
Clone_based_ensembl_gene
Clone_based_ensembl_transcript
RFAM_trans_name
miRBase_trans_name
EntrezGene_trans_name);

  push @list, $dbname."_trans_name";
  push @list, $dbname;

  my $sth = $dbi->prepare("select source_id from source where name like ?");
  
  my $source_error = 0;
  foreach my $source (@list){
    my $id  = undef;
    $sth->execute($source);
    $sth->bind_columns(\$id);
    $sth->fetch();
    if(!defined($id)){
      warn "Could not find external database name $source\n";
      $source_error++;
    }
    else{
      $dbname_to_source_id{$source} = $id;
    }
  }
  if($source_error){
    carp "Could not find name for $source_error database name.\nTherefore Exiting.\nPlease add these sources";
  }
  return \%dbname_to_source_id;
}

sub delete_old_data{
  my ($self, $dbname_to_source_id, $dbi) = @_;

  my $dbname = $self->get_official_name();

 my @sources = qw(
Clone_based_ensembl_gene
Clone_based_ensembl_transcript
EntrezGene_trans_name
RFAM_trans_name
miRBase_trans_name);

  push @sources, $dbname."_trans_name";  

  my @source_ids = map {$dbname_to_source_id->{$_}} @sources;
  my $list = join(", ",@source_ids);
  

  print "LIST to delete $list\n";
					    
       
  my $sql =(<<"DE1");
DELETE s
  FROM synonym s, xref x 
    WHERE s.xref_id = x.xref_id AND 
          x.source_id in ( $list );
DE1

  my $sth = $dbi->prepare($sql);
  $sth->execute();

 
  my $del_identity_sql =(<<"DE2");
DELETE i 
  FROM object_xref o, xref x, identity_xref i
    WHERE i.object_xref_id = o.object_xref_id AND
           x.xref_id = o.xref_id AND
            x.source_id in ( $list )
DE2
  $sth = $dbi->prepare($del_identity_sql);
  $sth->execute();
 
  my $del_ox_sql = (<<"DE3");
DELETE o 
  FROM object_xref o, xref x 
    WHERE x.xref_id = o.xref_id AND
           x.source_id in ( $list )
DE3
  $sth = $dbi->prepare($del_ox_sql);
  $sth->execute();
 
  my $del_x_sql = "delete x from xref x where x.source_id in ( $list )";

  $sth = $dbi->prepare($del_x_sql);
  $sth->execute();
  return;
}


sub reset_display_xrefs{
  my $self = shift;
  my $dbi = shift;

  my $sth =  $dbi->prepare("update transcript_stable_id set display_xref_id = null");
  $sth->execute;

  $sth = $self->xref->dbc->prepare("UPDATE gene_stable_id SET display_xref_id = null, desc_set =0");
  $sth->execute;

  return;
}

1;

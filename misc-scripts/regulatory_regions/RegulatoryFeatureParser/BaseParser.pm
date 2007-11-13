package RegulatoryFeatureParser::BaseParser;

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Analysis;

# Base functionality for external_feature parsers

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = {};

  #my $self = {(
#			  feature_sets => {(
#								vista   => {('VISTA enhancer set' => undef)},
#								cisred  => {(
#											 'cisRED search regions' => undef,
#											 'cisRED group motifs'   => undef,
#											)},
#								miranda => {('miRanda miRNA' => undef)},
#							   )},
#			  )};
  bless $self, $class;

  #validate and set type, analysis and feature_set here
  my ($type, $db, $clobber, $archive) = rearrange(['TYPE', 'DB', 'CLOBBER', 'ARCHIVE'], @_);
  
  throw('You must define a type of external_feature to import') if(! defined $type);

  if (! ($db && ref($db) &&
		 $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'))){
	throw('You must provide a valid Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  }

  throw('You can only specify either -clobber or -archive, but not both') if($clobber && $archive);


  $self->{'db'} = $db;
  $self->{type} = $type;
  $self->{'clobber'} = $clobber if defined $clobber;
  $self->{'archive'} = $archive if defined $archive;
  
  print ":: Parsing and loading $type ExternalFeatures";

  return $self;

}

sub db{
  my ($self, $db) = @_;

  if($db){

	if(! $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor')){
	  throw('You must prove a Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
	}
	
	$self->{'db'} = $db;
  }

  return $self->{'db'};
}


sub set_feature_sets{
  my $self = shift;

  throw('Must provide a set feature_set config hash') if ! defined $self->{'feature_sets'};


  my $fset_adaptor = $self->db->get_FeatureSetAdaptor;
  my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
  
  foreach my $fset_name(keys %{$self->{feature_sets}}){
	
	#check feature_type is validated?

	my $fset = $fset_adaptor->fetch_by_name($fset_name);

	#what about data sets?
	#we don't need data sets for external_sets!

	if(defined $fset){

	  if($self->{'clobber'}){
		#Need to clobber any DBEntries first!!!

		if($self->{'type'} eq 'cisred'){
		  my @ext_feat_ids = 	@{$self->db->dbc->db_handle->selectall_arrayref('select external_feature_id from external_feature where feature_set_id='.$fset->dbID)};
		  
		  if(@ext_feat_ids){

			my ($core_ext_dbid) = @{$self->db->dbc->db_handle->selectall_arrayref('select external_db_id from external_db where db_name="core"')};

			if($core_ext_dbid){
			  #double table delete?
			  my $sql = "delete x, ox from object_xref ox, xref x where ox.ensembl_object_type='ExternalFeature' and x.external_db_id=$core_ext_dbid and ox.xref_id=x.xref_id and ensembl_id in(".join(', ', @ext_feat_ids).')';
			  print ":: Clobbering xrefs for $fset_name\n";
			  $self->db->dbc->do($sql);
			}
		  }
	
		  print ":: Clobbering old features for external feature_set:\t$fset_name\n";
		  my $sql = 'delete from external_feature where feature_set_id='.$fset->dbID;
		  $self->db->dbc->do($sql);
		}

	  }
	  elsif($self->{'archive'}){
		my $archive_fset =  $fset_adaptor->fetch_by_name($fset_name."_v".$self->{'archive'});

		if(defined $archive_fset){
		  throw("You are trying to create an archive external feature_set which already exists:\t${fset_name}_v".$self->{archive});
		}

		my $sql = "UPDATE feature_set set name='$fset_name}_v".$self->{archive}."' where name='$fset_name'";
		$self->db->dbc->do($sql);
		undef $fset;
	  }else{
		throw("You are trying to create an external feature_set which already exists:\t$fset_name\nMaybe to want to clobber or archive?");
	  }
	}

	if(!defined $fset){
	  #don't need to use RNAFeatureType here as this is the setwide generic feature_type
	  #or do we have separate tables for external_feature and external_rna_feature?

	  #validate analysis first
	  my $analysis = $analysis_adaptor->fetch_by_logic_name($self->{'feature_sets'}{$fset_name}{'analysis'}{'-logic_name'});
	  if(! defined $analysis){
		
		print ':: Analysis '.$self->{'feature_sets'}{$fset_name}{'analysis'}{'-logic_name'}.
		  " not found, storing from config hash\n";
		
		$analysis_adaptor->store(Bio::EnsEMBL::Analysis->new
								 (
								  %{$self->{'feature_sets'}{$fset_name}{'analysis'}}
								  #-logic_name    => $self->{'feature_sets'}{$fset_name}{'analysis'}{'logic_name'},
								  #-description   => $self->{'feature_sets'}{$fset_name}{'analysis'}{'description'},
								  #-display_label => $self->{'feature_sets'}{$fset_name}{'analysis'}{'display_label'},
								  #-diplayable    => $self->{'feature_sets'}{$fset_name}{'analysis'}{'displayable'},
								 ));
		
		$analysis = $analysis_adaptor->fetch_by_logic_name($self->{'feature_sets'}{$fset_name}{'analysis'});
	  }

	  #replace hash config with object
	  $self->{'feature_sets'}{$fset_name}{'analysis'} = $analysis;

	  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
													 -name         => $fset_name,
													 -type         => 'external',
													 -analysis     => $self->{'feature_sets'}{$fset_name}{'analysis'},
													 -feature_type => ${$self->{'feature_sets'}{$fset_name}{'feature_type'}},
													);

	  ($fset) = @{$self->db->get_FeatureSetAdaptor->store($fset)};
	}

	#Now replace config hash with object
	$self->{feature_sets}{$fset_name} = $fset;
  }

  return;
}



# --------------------------------------------------------------------------------
# Delete existing regulatory features etc

#This needn't be done in line any more as we can just make the old set not displayable

sub delete_existing {

  my ($db_adaptor, $type) = @_;


  throw('delete existing is deprecated');

  my $t = lc($type);

  print "Deleting existing features & related data for type $type\n";

  # Delete any regulatory_feature_coding entries first
  my $sth = $db_adaptor->dbc->prepare("DELETE rft FROM regulatory_feature rfeat, regulatory_factor_coding rft, analysis a WHERE rfeat.regulatory_factor_id=rft.regulatory_factor_id AND a.analysis_id=rfeat.analysis_id AND LOWER(a.logic_name)=?");
  $sth->execute($t);

  # now delete interlinked regulatory_feature, regulatory_factor and regulatory_feature_object entries
  $sth = $db_adaptor->dbc->prepare("DELETE rfeat, rfact, rfo FROM regulatory_feature rfeat, regulatory_factor rfact, regulatory_feature_object rfo, analysis a WHERE rfeat.regulatory_feature_id=rfo.regulatory_feature_id AND rfeat.regulatory_factor_id=rfact.regulatory_factor_id AND rfeat.analysis_id=a.analysis_id AND LOWER(a.logic_name)=?");
  $sth->execute($t);

  # delete dangling regulatory_factors
  $sth = $db_adaptor->dbc->prepare("DELETE rfact FROM regulatory_feature rfeat, regulatory_factor rfact, analysis a WHERE rfeat.regulatory_factor_id=rfact.regulatory_factor_id AND a.analysis_id=rfeat.analysis_id AND LOWER(a.logic_name)=?");
  $sth->execute($t);

  # and finally any dangling regulatory features
  $sth = $db_adaptor->dbc->prepare("DELETE rfeat FROM regulatory_feature rfeat, analysis a WHERE a.analysis_id=rfeat.analysis_id AND LOWER(a.logic_name)=?");
  $sth->execute($t);

  # Delete search regions; they have a different analysis_id
  if ($type eq "cisred") {
    my $sr_type = $type . "_search";
    die "Can't find analysis for $sr_type " unless validate_type($db_adaptor, $sr_type);
    my $anal_sth = $db_adaptor->dbc->prepare("SELECT analysis_id FROM analysis WHERE LOWER(logic_name)=?");
    $anal_sth->execute($sr_type);
    my $anal = ($anal_sth->fetchrow_array())[0];
    $sth = $db_adaptor->dbc->prepare("DELETE FROM regulatory_search_region WHERE analysis_id=?");
    $sth->execute($anal);
  }

}

# --------------------------------------------------------------------------------
# Check that the type specified corresponds to an analysis

sub validate_and_store_feature_types{
  my $self = shift;

  # LBL enhancers have both positive and negative varieties

  #feature type class/naming is logical but not intuitive
  #+-----------------+-------------------------+----------+----------------------------------------------------+
  #| feature_type_id | name                    | class    | description                                        |
  #+-----------------+-------------------------+----------+----------------------------------------------------+
  #|          398680 | VISTA Enhancer          | Enhancer | Enhancer identified by positive VISTA assay        |
  #|          398681 | VISTA Target - Negative | Region   | Enhancer negative region identified by VISTA assay |
  #+-----------------+-------------------------+----------+----------------------------------------------------+


  #if (lc($type) eq 'VISTA') {
  #  return (validate_type($db_adaptor, 'VISTA Enhancer') && validate_type($db_adaptor, 'VISTA Target - Negative'));
  #}

  #my $sth = $self->db->dbc->prepare("SELECT analysis_id FROM analysis WHERE logic_name=?");

  

  #remove lc as mysql doesn't care about case
  #$sth->execute($type);
  #if ($sth->fetchrow_array()) {
  #  print "Type $type is valid\n";
  #} else {
  #  print "Type $type is not valid - is there an entry for $type in the analysis table?\n";
  #  return 0;
  #}

  my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;

  foreach my $ftype_name(keys %{$self->{'feature_types'}}){

	my $ftype = $ftype_adaptor->fetch_by_name($ftype_name);


	if(! defined $ftype){
	  print ":: FeatureType '".$ftype_name."' for external feature_set ".$self->{'type'}." not present\n".
		":: Storing using type hash definitions\n";
	
	  $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new(
													   -name => $ftype_name,
													   -class => $self->{'feature_types'}{$ftype_name}{'class'},
													   -description => $self->{'feature_types'}{$ftype_name}{'description'},
													  );
	  ($ftype) = @{$ftype_adaptor->store($ftype)};
	}

	#Replace hash config with object
	$self->{'feature_types'}{$ftype_name} = $ftype;
  }

  return;
}


# --------------------------------------------------------------------------------
# Find the maximum existing ID in a table.
sub find_max_id {

  my ($self, $table) = @_;

  my $row = @{$self->dbc->db_handle->selectall_arrayref("SELECT MAX(${table}_id) FROM $table")}[0];
  my $max_id = @{$row}[0];

  if (!defined $max_id) {
    print "No existing ${table}_ids, will start from 1\n";
    $max_id = 0;
  } else {
    print "Maximum existing ${table}_id = $max_id\n";
  }

  return $max_id;

}

# --------------------------------------------------------------------------------
# Build a cache of ensembl stable ID -> internal ID
# Return hashref keyed on {$type}{$stable_id}
# Type is always all lower case

sub build_display_name_cache {

  my ($self, @types) = @_;

  my %display_name_cache;

  #validate types here?

  foreach my $type (@types){#'gene', 'transcript', 'translation') { # Add exon here if required

	#No display_xref for translations!!
	throw('Cannot build cache for tranlsations, just use stable_id, check your implementatio') if $type eq 'translation';

    print ":: Caching stable ID ->  display_name links for ${type}s\n";

    my $count = 0;
	my $sql;

	print "dnadb is ".Data::Dumper::Dumper($self->db->dnadb->dbc);
	
	#could do a fetchall_hashref here?
	my $sth= $self->db->dnadb->dbc->prepare("SELECT t.stable_id, x.display_label FROM s.${type}_stable_id, t.${type} left join x.xref on t.display_xref_id=x.xref_id and s.${type}_id=t.${type}_id");

	

    $sth->execute();
    my ($stable_id, $display_name);
    $sth->bind_columns(\$stable_id, \$display_name);

    while ($sth->fetch) {
	  #removed type
      $display_name_cache{$stable_id} = $display_name;
      $count++ if $display_name;#will this be ! NULL?
    }

    print ":: Got $count $type stable ID -> display_name mappings\n";

  }

  return %display_name_cache;
}

# --------------------------------------------------------------------------------
# Upload features and factors etc to database

sub upload_features_and_factors {

  my ($self, $objects) = @_;


  throw("upload_features_and_factors is deprecated");

  my $dbc = $self->db->dbc;



  # we need to convert this to use the external_feature API
  #this means we need to build external_features in the other modules and test for each RNA_feature_type/motif_feature_type?
  #keep as feature_types for now, maybe next release?

  my $feature_sth = $dbc->prepare("INSERT INTO regulatory_feature (regulatory_feature_id, name, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id, regulatory_factor_id) VALUES(?,?,?,?,?,?,?,?)");
  my $factor_sth = $dbc->prepare("INSERT INTO regulatory_factor (regulatory_factor_id, name, type) VALUES(?,?,?)");
  my $feature_object_sth = $dbc->prepare("INSERT INTO regulatory_feature_object (regulatory_feature_id, ensembl_object_type, ensembl_object_id, influence, evidence) VALUES(?,?,?,?,?)");

  my $sr_sth = $dbc->prepare("INSERT INTO regulatory_search_region (name, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, ensembl_object_type, ensembl_object_id, analysis_id ) VALUES(?,?,?,?,?,?,?,?)");

  print "Uploading " . scalar(@{$objects->{FEATURES}}) . " features ...\n";

  foreach my $feature (@{$objects->{FEATURES}}) {

    $feature_sth->execute($feature->{INTERNAL_ID},
			  $feature->{NAME},
			  $feature->{SEQ_REGION_ID},
			  $feature->{START},
			  $feature->{END},
			  $feature->{STRAND},
			  $feature->{ANALYSIS_ID},
			  $feature->{FACTOR_ID});

    if ($feature->{ENSEMBL_TYPE} && $feature->{ENSEMBL_ID}) {
      $feature_object_sth->execute($feature->{INTERNAL_ID},
				   $feature->{ENSEMBL_TYPE},
				   $feature->{ENSEMBL_ID},
				   $feature->{INFLUENCE},
				   $feature->{EVIDENCE});
    }

  }

  if ($objects->{FACTORS}) {

    print "Uploading " . scalar(@{$objects->{FACTORS}}) . " factors ...\n";

    foreach my $factor (@{$objects->{FACTORS}}) {

      $factor_sth->execute($factor->{INTERNAL_ID},
			   $factor->{NAME},
			   $factor->{TYPE});

    }

  }

  if ($objects->{SEARCH_REGIONS}) {

    print "Uploading " . scalar(@{$objects->{SEARCH_REGIONS}}) . " search regions ...\n";

    foreach my $search_region (@{$objects->{SEARCH_REGIONS}}) {

      $sr_sth->execute($search_region->{NAME},
		       $search_region->{SEQ_REGION_ID},
		       $search_region->{START},
		       $search_region->{END},
		       $search_region->{STRAND},
		       $search_region->{ENSEMBL_OBJECT_TYPE},
		       $search_region->{ENSEMBL_OBJECT_ID},
		       $search_region->{ANALYSIS_ID});


    }

  }

  print "Done\n";

}

# --------------------------------------------------------------------------------
# Project a feature from one slice to another
sub project_feature {
  my ($self, $feat, $new_assembly) = @_;

  # just use a SimpleFeature for convenience
  #my $feat = Bio::EnsEMBL::SimpleFeature->new
  #  (-start    => $start,
  #   -end      => $end,
  #   -strand   => $strand,
  #   -slice    => $slice,
  #   -analysis => $analysis,
  #   -display_label => $label,
  #   -score   => 0);

  # project feature to new assembly
  my $feat_slice = $feat->feature_Slice;
  my @segments = @{ $feat_slice->project('chromosome', $new_assembly) };


  #next what?
  #next unless (@segments);
  #next if (scalar(@segments) > 1);

  if(! @segments){
	print "Failed to project feature:\t".$feat->display_label."\n";
  }
  elsif(scalar(@segments) >1){
	print "Failed to project feature to distinct location:\t".$feat->display_label."\n";
	return;
  }

  my $proj_slice = $segments[0]->to_Slice;

  if($feat_slice->length != $proj_slice->length){
	print "Failed to project feature to comparable length region:\t".$feat->display_label."\n";
	return;
  }


  # everything looks fine, so adjust the coords of the feature
  $feat->start($proj_slice->start);
  $feat->end($proj_slice->end);
  $feat->strand($proj_slice->strand);
  my $slice_new_asm = $self->slice_adaptor->fetch_by_region('chromosome', $proj_slice->seq_region_name, undef, undef, undef, $new_assembly);
  $feat->slice($slice_new_asm);

  return $feat;

}

sub slice_adaptor{
  my $self = shift;

  if(! defined $self->{'slice_adaptor'}){
	$self->{'slice_adaptor'} = $self->db->get_SliceAdaptor;
  }
  
  return $self->{'slice_adaptor'};
}


# --------------------------------------------------------------------------------

#sub get_blank_factor_id () {

#  my ($self, $db_adaptor) = @_;

#  my $sth = $db_adaptor->dbc->prepare("SELECT regulatory_factor_id FROM regulatory_factor WHERE name=''");
#  $sth->execute();

#  my ($factor_id) = $sth->fetchrow_array();

#  if ($factor_id) {
#    print "Found existing blank factor, id = $factor_id\n";
#  } else {
#     $db_adaptor->dbc->do("INSERT INTO regulatory_factor (name) VALUES ('')");
#     $sth->execute();
#     ($factor_id) = $sth->fetchrow_array();
#     print "Created new blank factor, id = $factor_id\n";
#  }

#  return $factor_id;

#}


1;

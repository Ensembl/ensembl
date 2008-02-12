#Nath
#Implemented strand check and unmapped object for anti-sense mapping
#Changed transcript_probeset_count cache to use dbID instead of complete name
#Renamed check_names_match to validate_arrays and extended to perform probeset warning, external_db name guess and returns cache, user array param validation and association check
#Merged delete_unampped objects into delete_existing_xrefs due to dependency. Now throws if arrays param set until force_delete specified, deletes using IN list
#cache_arrays_per_probeset now uses arrays param to reduce size of cache
#removed unmapped_object dbID = 2000>
#Change all probe specific unmapped object to reeflect the individual probe rather than the probeset
#Updated logs
#Updated docs
#Added control of promiscuous probesets and unmapped objects


use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;


#$| = 1; # auto flush stdout

my ($transcript_host, $transcript_user, $transcript_pass, $transcript_dbname,
    $oligo_host, $oligo_user, $oligo_pass, $oligo_dbname, $five_utr, $three_utr,
    $xref_host, $xref_user, $xref_pass, $xref_dbname, $force_delete,
    $max_transcripts, @arrays, $delete, $no_triage, $health_check);

my ($oligo_db, $xref_db, %promiscuous_probesets, %transcripts_per_probeset, @unmapped_objects, $um_obj,
	%transcript_ids , %transcript_probeset_count, %arrays_per_probeset, %array_probeset_sizes);

# Default options
my $transcript_port = 3306; 
my $oligo_port = 3306; 
my $xref_port = 3306;

my $max_mismatches = 1;
my $utr_length = 2000;
my $max_transcripts_per_probeset = 100;
my $mapping_threshold = 0.5;

GetOptions(
		   'transcript_host=s'      => \$transcript_host,
           'transcript_user=s'      => \$transcript_user,
           'transcript_port=i'      => \$transcript_port,
           'transcript_pass=s'      => \$transcript_pass,
           'transcript_dbname=s'    => \$transcript_dbname,
		   'oligo_host=s'           => \$oligo_host,
           'oligo_user=s'           => \$oligo_user,
           'oligo_port=i'           => \$oligo_port,
           'oligo_pass=s'           => \$oligo_pass,
           'oligo_dbname=s'         => \$oligo_dbname,
		   'xref_host=s'            => \$xref_host,
           'xref_user=s'            => \$xref_user,
           'xref_port=i'            => \$xref_port,
           'xref_pass=s'            => \$xref_pass,
           'xref_dbname=s'          => \$xref_dbname,
		   'mismatches=i'           => \$max_mismatches,
           'utr_length=s'           => \$utr_length,
		   'max_probesets=i'        => \$max_transcripts_per_probeset,
		   'max_transcripts=i'      => \$max_transcripts,
		   'threshold=s'            => \$mapping_threshold,
		   'arrays=s'               => \@arrays,
		   'delete'                 => \$delete,
		   'force_delete'           => \$force_delete,
		   'no_triage'              => \$no_triage,
		   'health_check'           => \$health_check,
           'help'                   => sub { usage(); exit(0); }
		  );

@arrays = split(/,/,join(',',@arrays));#?

if(($utr_length =~ /\D/) && ($utr_length ne 'annotated')){
  die("Invalid utr_length parameter($utr_length).  Must be a number or 'annotated'");
}
else{
  $three_utr = $utr_length;
}

usage() if(!$transcript_user || !$transcript_dbname || !$transcript_host);


print 'Running on probe2trascript.pl on: '.`hostname`."\n";

my $transcript_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $transcript_host,
						       '-port'   => $transcript_port,
						       '-user'   => $transcript_user,
						       '-pass'   => $transcript_pass,
						       '-dbname' => $transcript_dbname);

if ($oligo_host && $oligo_dbname && $oligo_user) {

  $oligo_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $oligo_host,
						 '-port'   => $oligo_port,
						 '-user'   => $oligo_user,
						 '-pass'   => $oligo_pass,
						 '-dbname' => $oligo_dbname);

} else {

  print "No oligo database specified, defaulting to $transcript_host:$transcript_port:$transcript_dbname\n";
  $oligo_db = $transcript_db;

}


if ($xref_host && $xref_dbname && $xref_user) {

  $xref_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $xref_host,
						'-port'   => $xref_port,
						'-user'   => $xref_user,
						'-pass'   => $xref_pass,
						'-dbname' => $xref_dbname);

} else {

  print "No xref database specified, writing xrefs to $transcript_host:$transcript_port:$transcript_dbname\n";
  $xref_db = $transcript_db;

}

#This validates arrays
#Why would we ever want to write the xrefs to a different DB?
my %array_name_cache =  %{&validate_arrays($oligo_db, $xref_db)};
$delete = $force_delete if $force_delete;

#Merge these as they are related and we need to force unmapped check first

delete_existing_xrefs($xref_db) if ($delete);
check_existing_and_exit($oligo_db, $xref_db);#???

#Do we need this healthcheck mode now?  We're doing it all anyway and fails if any tests fail
if ($health_check){
  print "Healthcheck passed\n";
  exit 0;
}

my $transcript_adaptor = $transcript_db->get_TranscriptAdaptor();
my $slice_adaptor = $transcript_db->get_SliceAdaptor();
my $oligo_feature_adaptor = $oligo_db->get_OligoFeatureAdaptor();
my $db_entry_adaptor = $xref_db->get_DBEntryAdaptor();
my $analysis_adaptor = $xref_db->get_AnalysisAdaptor();
my $unmapped_object_adaptor = $xref_db->get_UnmappedObjectAdaptor();

my $analysis = get_or_create_analysis($analysis_adaptor);

my $i = 0;
my $last_pc = -1;

my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();

my @transcripts = @{$transcript_adaptor->fetch_all()};
my $total =  scalar(@transcripts);

open (LOG, ">${transcript_dbname}_probe2transcript.log");

print "Identified ".scalar(@transcripts)." transcripts for probe mappinng\n";
print "Mapping, percentage complete: ";

foreach my $transcript (@transcripts) {

  my $pc = int ((100 * $i) / $total);

  if ($pc > $last_pc) {
    print "$pc ";
    $last_pc = $pc;
  }

  $i++;

  last if ($max_transcripts && $i >= $max_transcripts);

  my $stable_id = $transcript->stable_id();
  $transcript_ids{$stable_id} = $transcript->dbID(); # needed later

  if($utr_length eq 'annotated'){
	#$five_utr = Do we need to implement this?
	my $utr = $transcript->three_prime_utr;
	$three_utr = (defined $utr) ? $utr->length : 2000;

	#Should we set UTR to 0 if no UTR?
	#Should we default to 2000 here if UTR is less?
	#Probably dependent on quality of gene build
  }

  my $slice = $transcript->feature_Slice(); # note not ->slice() as this gets whole chromosome!
  #my $extended_slice = $slice->expand(0, $utr_length); # this takes account of strand
  my $extended_slice = $slice->expand(0, $three_utr); # this takes account of strand
  

  my $exons = $transcript->get_all_Exons();

  $rr->flush();


  foreach my $exon (@$exons) {
    my $start = $exon->start();
    my $end = $exon->end();
    $rr->check_and_register("exonic", $start, $end);
  }

  if ($transcript->strand() == 1) {
    $rr->check_and_register("utr", $extended_slice->end()-$three_utr, $extended_slice->end());
  } else {
    $rr->check_and_register("utr", $extended_slice->start(), $extended_slice->start() + $three_utr);
  }

  my $oligo_features;

 
  if (@arrays) {
    $oligo_features = $extended_slice->get_all_OligoFeatures(@arrays);
  } else {
    $oligo_features = $extended_slice->get_all_OligoFeatures();
  }

 
  #This works on the assumption that probesets are identical between arrays
  #i.e. if a probe set is present on different arrays, their probes are identical.
 
  foreach my $feature (@$oligo_features) {
  	#Here we need to skip the assignment if we have already seen the dbID for this transcript
	#Actually we need to count each dbID mapping
	#Then if we get some which don't match we need to resolve the differences 
	#between the array mappings and create 2 xrefs
	#Then change how arrays per probeset works
	#What about if we have > 1 mapping for a given probe?
	#Can we use this in the xref desc?


	my $probe = $feature->probe();
    #my $probe_name = ${$probe->get_all_complete_names()}[0];
	#No need to use this as a unique key as we can use dbID, this wouldn't work anyway without redundant dbID records

	#my $probe_name = $probe->get_probename;
	my $dbID = $probe->dbID;
	my $probeset = $probe->probeset();


	#We should only get one probe for all arrays
	#So if we get two we know we have >1 mapping
	#So we need to test again.

	if ($transcript->strand() != $feature->strand()){
	  print LOG "Unmapped anti-sense " . $transcript->stable_id . "\t${probeset}\tdbID:${dbID}\n";
	  
	  if (!$no_triage) {
		
		#Use of internal dbID in identifier is dodge
		#can we use the actual probe names here?
		$um_obj = new Bio::EnsEMBL::UnmappedObject(-type       => 'probe2transcript',
													  -analysis   => $analysis,
													  -identifier => "${probeset}:${dbID}",
													  -summary    => "Unmapped anti-sense",
													  -full_desc  => "Probe mapped to opposite strand of transcript",
													  -ensembl_object_type => 'Transcript',
													  -ensembl_id => $transcript->dbID());
		
		&cache_and_load_unmapped_objects($um_obj);

		next;
	  }
	}


	#Are we dealing with each probeset for each array individually?
	#Can we generalise across arrays?
	#We need to deal with the probe dbIDs in the cache?
	#Are the probe names being used for anything?
	#We're are we deconvoluting the probe to array relationship?

    
    my $transcript_key = $transcript->stable_id() . ":" . $probeset;
    my $probe_length = $probe->probelength();
    my $min_overlap = ($probe_length - $max_mismatches);
    my $exon_overlap = $rr->overlap_size("exonic", $feature->seq_region_start(), $feature->seq_region_end());
    my $utr_overlap  = $rr->overlap_size("utr",    $feature->seq_region_start(), $feature->seq_region_end());

    if ($exon_overlap >= $min_overlap) {
	  #	  print "Exon overlap $exon_overlap excedes minimum overlap $min_overlap\n";
      $transcript_probeset_count{$transcript_key}{$dbID}++;
    } 
	elsif ($utr_overlap >= $min_overlap) {
	  #	  print "UTR overlap $utr_overlap excedes minimum overlap $min_overlap\n";
      $transcript_probeset_count{$transcript_key}{$dbID}++;
    } 
	else { # must be intronic
	  #This is for an individual probe!
	  print LOG "Unmapped intronic " . $transcript->stable_id . "\t${probeset}\tdbID:${dbID}\n";
	  
	  if (!$no_triage) {
		
		$um_obj = new Bio::EnsEMBL::UnmappedObject(-type       => 'probe2transcript',
													  -analysis   => $analysis,
													  -identifier => "${probeset}:${dbID}",
													  -summary    => "Unmapped intronic",
													  -full_desc  => "Probe mapped to intronic region of transcript",
													  -ensembl_object_type => 'Transcript',
													  -ensembl_id => $transcript->dbID());
		
		&cache_and_load_unmapped_objects($um_obj);

	
	  }
	}
  }

 
  # TODO - make external_db array names == array names in OligoArray!
  # change oligo_array.name to be oligo_array.external_db_id ?

}

print "\n";

cache_arrays_per_probeset($oligo_db);

print "Writing xrefs\n";
my $um_cnt = 0;

# now loop over all the mappings and add xrefs for those that have a suitable number of matches
foreach my $key (keys %transcript_probeset_count) {

  my ($transcript, $probeset) = split (/:/, $key);

  # store one xref/object_xref for each array-probeset-transcript combination
  foreach my $array (split (/ /, $arrays_per_probeset{$probeset})) {	
	#Not sure about this?
	#If we're restricting the an array list, then we're ignoring the genuine array relationship
	#Which may give spurious results
	#We need to think how we're rolling back xrefs and unmapped objects anyway
	#As these are done one a probe_set/probe level, not an array level
	#So we mayve have duplicates if we don't clean
	#But we have no way of cleaning on an array level without deleting all
	#Can we clean on a probe/probeset level?


	#This should never happen now, as we're enforcing running on all associated arrays
    next if (@arrays && find_in_list($array, @arrays,) == -1 );

    my $size_key = $array . ":" . $probeset;
    my $probeset_size = $array_probeset_sizes{$size_key};


	#Is this the only place we're using the complete name?
    my $hits = scalar(keys %{$transcript_probeset_count{$key}});

    if ($hits / $probeset_size >= $mapping_threshold) {
	  #This is inc'ing an undef?
	  $transcripts_per_probeset{$probeset}++;

      # only create xrefs for non-promiscuous probesets

      #XXX
      #add_xref($transcript_ids{$transcript}, $probeset, $db_entry_adaptor, $array_name_cache{$array}, $probeset_size, $hits);
      #print LOG "$probeset\t$transcript\tmapped\t$probeset_size\t$hits\n";

	  if ($transcripts_per_probeset{$probeset} <= $max_transcripts_per_probeset) {     
        add_xref($transcript_ids{$transcript}, $probeset, $db_entry_adaptor, $array, $probeset_size, $hits);
        print LOG "$probeset\t$transcript\tmapped\t$probeset_size\t$hits\n";
     
      } 
	  else {
        print LOG "$probeset\t$transcript\tpromiscuous\t$probeset_size\t$hits\tCurrentTranscripts".$transcripts_per_probeset{$probeset}."\n";
        push @{$promiscuous_probesets{$probeset}}, $transcript_ids{$transcript};
	  }

    } 
	else {

      print LOG "$probeset\t$transcript\tinsufficient\t$probeset_size\t$hits\n";
	  
	  if (!$no_triage) {

		#Can/should we concentrate all unmapped info into one record
		#Currently getting one for each probe and each probeset
		
		$um_obj = new Bio::EnsEMBL::UnmappedObject(-type       => 'probe2transcript',
												   -analysis   => $analysis,
												   -identifier => $probeset,
												   -summary    => "Insufficient hits",
												   -full_desc  => "Probeset had an insufficient number of hits (probeset size = $probeset_size, hits = $hits)",
												   -ensembl_object_type => 'Transcript',
												   -ensembl_id => $transcript_ids{$transcript});
		

		
		&cache_and_load_unmapped_objects($um_obj);
      }
    }
  }
}



#Now update promiscuous probesets
print "Updating ".scalar(keys %promiscuous_probesets)."\n";

foreach my $probeset(keys %promiscuous_probesets){

  #First delete mapped object_xrefs
  #As there is a chance that probes might be xreffed to a non-transcript entity
  #Deleting ox and x at the same time would orphan any non-transcript ox's
  $xref_db->dbc()->do("DELETE ox FROM object_xref ox, xref x WHERE x.dbprimary_acc='$probeset' AND x.xref_id=ox.xref_id AND ensembl_object_type='Transcript'");

  #Any other oxs?
  if(! @{ $xref_db->dbc->db_handle->selectall_arrayref("SELECT ox.object_xref_id from object_xref ox, xref x WHERE x.dbprimary_acc='$probeset' AND x.xref_id=ox.xref_id")}){
	#Then delete xref
	 $xref_db->dbc()->do("DELETE FROM xref WHERE dbprimary_acc='$probeset'");
  }
  

  #Now load all unmapped objects
  #One for all arrays rather than one for each
  foreach my $transcript(@{$promiscuous_probesets{$probeset}}){

	$um_obj = new Bio::EnsEMBL::UnmappedObject
	  (
	   -type       => 'probe2transcript',
	   -analysis   => $analysis,
	   -identifier => $probeset,
	   -summary    => 'Promiscuous probeset',
	   -full_desc  => 'Probeset maps to '.$transcripts_per_probeset{$probeset}.' transcripts (max 100)',
	   -ensembl_object_type => 'Transcript',
	   -ensembl_id => $transcript_ids{$transcript}
	  );
		
	&cache_and_load_unmapped_objects($um_obj);
  }
}

# Find probesets that don't match any transcripts at all, write to log file
log_orphan_probes();

close (LOG);

# upload triage information if required
if (!$no_triage) {

  print "Uploading unmapped reasons to xref database\n";
  $unmapped_object_adaptor->store(@unmapped_objects);
  $um_cnt += scalar(@unmapped_objects);
  print "Loaded a total of $um_cnt UnmappedObjects to the DB\n";
}

# ----------------------------------------------------------------------

# only loads unless cache hits size limit

sub cache_and_load_unmapped_objects{
  my @um_obj = @_;

  push @unmapped_objects, @um_obj;

  if(scalar(@unmapped_objects) >10000){
	#This is setting the dbID of the unmapped obj, not the unmapped reason
	#Is this really working?
	#This won't do anything as the dbID attr is ignored on store
	#What was this trying to solve?
	#Was Martin mysqlimporting and overwriting data?
	#Surely this would fail on import with duplicate keys?
	#Remove for now
	#if($first_cache){
	#  $unmapped_objects[0]->dbID('2000');
	#  $first_cache = 0;
	#}

	$um_cnt += scalar(@unmapped_objects);

	$unmapped_object_adaptor->store(@unmapped_objects);
	@unmapped_objects = ();
  }
}





# ----------------------------------------------------------------------

sub log_orphan_probes {

  print "Logging probesets that don't map to any transcripts\n";

  foreach my $probeset (keys %arrays_per_probeset) {

    if (!$transcripts_per_probeset{$probeset}) {

      print LOG "$probeset\tNo transcript mappings\n";

      if (!$no_triage && $probeset) {
		$um_obj = new Bio::EnsEMBL::UnmappedObject(-type       => 'probe2transcript',
												   -analysis   => $analysis,
												   -identifier => $probeset,
												   -summary    => "No transcript mappings",
												   -full_desc  => "Probeset did not map to any transcripts");
		
		&cache_and_load_unmapped_objects($um_obj);
      }
    }
  }
}

# ----------------------------------------------------------------------

sub cache_arrays_per_probeset {
  my $db = shift;

  print "Caching arrays per probeset\n";
  my $sql;#do not need distinct on count as we're linking by array?

  if(@arrays){
	$sql = 'SELECT op.probeset, oa.name, count(op.oligo_probe_id) FROM oligo_probe op, oligo_array oa WHERE oa.oligo_array_id=op.oligo_array_id and oa.name in ("'.join('", "', @arrays).'") GROUP BY op.probeset, oa.name';
  }
  else{
	$sql = 'SELECT op.probeset, oa.name, count(op.oligo_probe_id) FROM oligo_probe op, oligo_array oa WHERE oa.oligo_array_id=op.oligo_array_id GROUP BY op.probeset, oa.name';
  }


  my $sth = $db->dbc()->prepare($sql);



  $sth->execute();
  my ($probeset, $array, $probeset_size);
  $sth->bind_columns(\$probeset, \$array, \$probeset_size);

  my $last_probeset = "";
  my $arrays = "";

  while($sth->fetch()){

    if ($probeset eq $last_probeset) {
      $arrays .= " " if ($arrays);
      $arrays .= $array;
    } else {
      $arrays = $array;
    }
    $arrays_per_probeset{$probeset} = $arrays;
    $last_probeset = $probeset;
    my $key = $array . ":" . $probeset;
    $array_probeset_sizes{$key} = $probeset_size;

  }

  $sth->finish();

}


# ----------------------------------------------------------------------

sub pc {

  my ($a, $total) = @_;

  return "?" if (!$total);

  my $number = 100 * $a / $total;
  my $pad = "";
  $pad .= " " if ($number < 100);
  $pad .= " " if ($number < 10);

  return $pad . sprintf "%3.2f", $number;

}

# ----------------------------------------------------------------------

sub add_xref {

  my ($transcript_id, $probeset, $dbea, $array, $probeset_size, $hits) = @_;

  my $txt = "probeset_size $probeset_size hits $hits";

  # TODO - this only works if external_db db_name == array name; currently needs
  # some manual hacking to acheive this

  my $dbe = new Bio::EnsEMBL::DBEntry( -adaptor              => $dbea,
				       -primary_id           => $probeset,
				       -version              => "1",
				       -dbname               => $array,
				       -release              => "1",
				       -display_id           => $probeset,
				       -description          => undef,
				       -primary_id_linkable  => 1,
				       -display_id_linkable  => 0,
				       -priority             => 1,
				       -db_display_name      => $array,
				       -info_type            => "MISC",  # TODO - change to PROBE when available
				       -info_text            => $txt);

  $dbea->store($dbe, $transcript_id, "Transcript");

}

# ----------------------------------------------------------------------

# Delete existing xrefs & object xrefs & unmapped objects. Use user-specified arrays if
# defined, otherwise all arrays.
# Now uses array name cache

sub delete_existing_xrefs {
  my $xref_db = shift;

  if(@arrays && ! $force_delete){
	die("You are attempting to delete all unampped objects even though you are only running a subsets of arrays.\n".
		  "This may result in losing unmapped information for other arrays.  If you really want to do this you must specify -force_delete")
  }

  #This deletes all unmapped objects, even if we're only (re-)running one array!
  print "Deleting ALL unmapped records\n";
  my $sql = 'DELETE a, ur, uo FROM analysis a, unmapped_reason ur, unmapped_object uo WHERE a.logic_name ="probe2transcript" AND a.analysis_id=uo.analysis_id AND uo.unmapped_reason_id=ur.unmapped_reason_id';
  $xref_db->dbc->do($sql);
	

  #Now delete xrefs for each array
  #can we change this to use an IN?
  #$del_sth = $dbentry_adaptor->dbc()->prepare("DELETE x, ox FROM xref x, object_xref ox, external_db e WHERE x.xref_id=ox.xref_id AND e.external_db_id=x.external_db_id AND e.db_name = ?");
  #foreach my $array (values %array_name_cache) {
  #  print "Deleting xrefs and object_xrefs for $array\n";
	#  $del_sth->execute($array);
  #}
  

  #Does this work? Will there no be some x's left?
  print "Deleting XREFs from:\t".join(' ', values %array_name_cache)."\n";
  $sql = 'DELETE x, ox FROM xref x, object_xref ox, external_db e WHERE x.xref_id=ox.xref_id AND e.external_db_id=x.external_db_id AND e.db_name in ("'.join('", "', values %array_name_cache).'")';
  $xref_db->dbc->do($sql);

  return;
}


# ----------------------------------------------------------------------

# Check if there are already xrefs defined, and exit if there are.
# Use user-specified arrays if defined, otherwise all arrays.
# Assumes external_db.dbname == oligo_array.name

sub check_existing_and_exit {

  my ($oligo_adaptor, $dbentry_adaptor) = @_;


  my $xref_sth = $dbentry_adaptor->dbc()->prepare("SELECT COUNT(*) FROM xref x, object_xref ox, external_db e WHERE x.xref_id=ox.xref_id AND e.external_db_id=x.external_db_id AND e.db_name = ?");

  foreach my $array (values %array_name_cache) {

    $xref_sth->execute($array);
    my @row2 = $xref_sth->fetchrow_array();
    if ($row2[0] > 0) {
      print "Array $array already has " . $row2[0] . " xrefs, exiting.\nThere may be other arrays with xrefs; use -delete to remove them if required.\n";
      exit(1);
    }

  }

  $xref_sth->finish();

}

# ----------------------------------------------------------------------

# Check that the array names in oligo_db appear in the external_db table in xref_db

#sub check_names_match {
sub validate_arrays{
 my ($oligo_adaptor, $dbentry_adaptor) = @_;

 my ($mismatch ,$sql);
 my %oligo_array_names;


 #Get all external_db entries and check probe_setsize
 #Name is not unique in core!!??
 $sql = "SELECT DISTINCT(name), probe_setsize FROM oligo_array";
 $sql .= ' where name in("'.join('", "', @arrays).'")' if @arrays;

 my $oligo_sth = $oligo_adaptor->dbc->prepare($sql);
 $oligo_sth->execute();

 print "Validating external_db names\n";

 while(my ($oa_name, $ps_size) = $oligo_sth->fetchrow_array()){

   my ($edb_name) = $dbentry_adaptor->dbc->db_handle->selectrow_array("SELECT db_name FROM external_db where db_name='${oa_name}'");

   #Try AFFY_ Can remove this when we change the naming convention
   if(! defined $edb_name){
	 ($edb_name) = $dbentry_adaptor->dbc->db_handle->selectrow_array("SELECT db_name FROM external_db where db_name='AFFY_${oa_name}'");
   }

   print "Array $oa_name does not have a probeset_size set, please rectify\n" if ! $ps_size;#boolean as we don't want 0 size


   #Should we add another try here using s/-/_/g?

   if( ! defined $edb_name){
	 print "Cannot find external_db for oligo_array:\t$oa_name\n";
	 $mismatch = 1;
   }
   else{
	 print "Found oligo_array\t$oa_name\tin external_db\t$edb_name\n";

	 $oligo_array_names{$oa_name} = $edb_name;
   }
 }

 $oligo_sth->finish();

 #Now do @arrays checks
 if(@arrays){
   print "Validating user specified arrays:\t@arrays\n";

   #Now check all specified arrays are valid
   if(scalar(@arrays) != scalar(keys %oligo_array_names)){

	 foreach my $array(@arrays){
	   
	   if(! exists $oligo_array_names{$array}){
		 print "Specified array does not exist in the oligo_array table:\t$array\n";
	   }
	 }
	 print "Could not find specified arrays in oligo_array table";
	 exit(1);
   }
   print "Found all external_dbs\n";

   #Now check we have all link arrays in @arrays
   #This cannot be guaranteed to work as some probesets may be missing on related arrays
   #Let's take a handfull here to reduce the risk of missing linked arrays
   # You souldn't really xref the arrays as subsets unless they are not related in anyway
   # In which case filtering @arrays in the cache will reduce the saize of the cache by removing unwanted array data
   # How do we tell the array groups?
   # 1 All old affy 3' and gneomics arrays should be mapped/xreffed together
   # 2 All new affy exon/gene ST arrays should be mapped/xreffed together
   # 3 Custom arrays can be mapped separately only if they have nor relation to pre-existing arays? Optional
   # Do we need another field in array? Design family? 

  
   #This is the ticket, will identify all linked arrays for a given array name
   my %linked_arrays;
   my $sth = $oligo_adaptor->dbc->prepare('select distinct(oa.name) from oligo_probe op, oligo_array oa where oa.oligo_array_id=op.oligo_array_id and op.oligo_probe_id in(select op.oligo_probe_id from oligo_probe op, oligo_array oa where oa.oligo_array_id=op.oligo_array_id and oa.name =?)');
   
   print "Getting all associated array data...this may take a while :|\n";
   
   foreach my $array(@arrays){

	 next if exists $linked_arrays{$array};

	 $sth->execute($array);

	 while(my ($larray) = $sth->fetchrow_array){
	   
	   $linked_arrays{$larray} = undef;
	 }
   }

   my @missing_arrays;

   foreach my $larray(keys %linked_arrays){

	 if(! grep(/^${larray}$/, @arrays)){
	   push @missing_arrays, $larray;
	 } 
   }

   if(@missing_arrays){
	 print "The array list you have specified is missing some associated arrays:\t@missing_arrays\n".
	   "This can cause incomplete xref data to be written.  Please check and start again\n";
	 exit(1);
   }
 }


 if ($mismatch) {
   print "At least one oligo array name was found that does not appear in external_db; this needs to be rectified before this script can be run successfully.\n";
   exit(1);
 }

 return \%oligo_array_names;

}



# ----------------------------------------------------------------------

# Find the index of an item in a list(ref), or -1 if it's not in the list.
# Only look for exact matches (case insensitive)

sub find_in_list {

  my ($item, @list) = @_;

  for (my $i = 0; $i < scalar(@list); $i++) {
    if (lc($list[$i]) eq lc($item)) {
      return $i;
    }
  }

  return -1;

}

# ----------------------------------------------------------------------

sub get_or_create_analysis {

  my ($analysis_adaptor) = @_;

  my $analysis = $analysis_adaptor->fetch_by_logic_name("probe2transcript");

  if (!$analysis) {

    my $id = $analysis_adaptor->store(new Bio::EnsEMBL::Analysis(-logic_name    => 'probe2transcript',
								 -program       => 'probe2transcript.pl',
								 -description   => 'Probe to transcript mapping',
								 -displayable   => '0'));

    $analysis = $analysis_adaptor->fetch_by_logic_name("probe2transcript");

  }

  return $analysis;

}

# ----------------------------------------------------------------------


sub usage {

  print << "EOF";

  Maps oligo probes to transcripts.

  perl $0 {options}

  Options ([..] indicates optional):

  READING TRANSCRIPTS:

   --transcript_host          The database server to read transcripts from.
		
   [--transcript_port]        The port to use for reading transcripts. Defaults to 3306.
		
   --transcript_user          Database username for reading transcripts.
		
   --transcript_pass          Password for transcript_user, if required.
		
   --transcript_dbname        Database name to read transcripts from.

  READING OLIGOS:

   --oligo_host          The database server to read oligo features from.
		
   [--oligo_port]        The port to use for reading oligo features. Defaults to 3306.
		
   --oligo_user          Database username for reading oligo featuress.
		
   --oligo_pass          Password for oligo_user, if required.
		
   --oligo_dbname        Database name to read oligo features from.

  WRITING XREFS:

   --xref_host          The database server to write xrefs to.
		
   [--xref_port]        The port to use for writing xrefs.. Defaults to 3306.
		
   --xref_user          Database username for xrefs. Must allow writing.
		
   --xref_pass          Password for xref_user, if required.
		
   --xref_dbname        Database name to write xrefs to.

  Note that if no oligo_host, xref_host etc is specified, oligo features will be read from,
  and xrefs written to, the database specified by the transcript_* parameters.

  Also, triage information will be written to the unmapped_object & unmapped_reason tables
  in the xref database, unless the -no_triage option is specified.

  GENERAL MAPPING OPTIONS:

  [--mismatches]      Allow up to this number of mismatches, inclusive.
                      Defaults to 1.
		
  [--utr_length]      Search this many bases downstream of the transcript
                      coding region as well. Defaults to 2000. Specify 'annotated' to use annotated lengths.
		
  [--max_probesets]   Don't store mappings to any 'promiscuous' probesets that map
                      to more than this number of transcripts. Defaults to 100.

  [--arrays]          Comma-separated list of arrays to use. Defaults to all arrays.

  [--threshold]       Fraction of probes per probeset that have to map. Default 0.5 

  MISCELLANEOUS:

  [--delete]          Delete existing xrefs and object_xrefs, and entries in unmapped_object.
                      No deletion is done by default.

  [--force_delete]    Forces deletion of all unmapped object info even if using a subset of arrays.

  [--max_transcripts] Only use this many transcripts. Useful for debugging.

  [--no_triage]       Don't write to the unmapped_object/unmapped_reason tables.

  [--health_check]    Only do sanity checks, then stop.

  [--help]            This text.


EOF

  exit(0);

}

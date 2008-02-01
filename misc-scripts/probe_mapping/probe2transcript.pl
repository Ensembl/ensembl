use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

my ($transcript_host, $transcript_port, $transcript_user, $transcript_pass, $transcript_dbname,
    $oligo_host, $oligo_port, $oligo_user, $oligo_pass, $oligo_dbname, $five_utr, $three_utr,
    $xref_host, $xref_port, $xref_user, $xref_pass, $xref_dbname,
    $max_mismatches, $utr_length, $max_transcripts_per_probeset, $max_transcripts, @arrays, $delete,
    $mapping_threshold, $no_triage, $health_check);

my $first_cache = 1;

GetOptions('transcript_host=s'      => \$transcript_host,
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
		   'no_triage'              => \$no_triage,
		   'health_check'           => \$health_check,
           'help'                   => sub { usage(); exit(0); });



# Default options
$transcript_port ||= 3306; $oligo_port ||= 3306; $xref_port ||= 3306;

$max_mismatches ||= 1;

$utr_length ||= 2000;

if(($utr_length =~ /\D/) && ($utr_length ne 'annotated')){
  die("Invalid utr_length parameter($utr_length).  Must be a number or 'annotated'");
}
else{
  $three_utr = $utr_length;
}


$max_transcripts_per_probeset ||= 100;

$mapping_threshold ||= 0.5;

@arrays = split(/,/,join(',',@arrays));

usage() if(!$transcript_user || !$transcript_dbname || !$transcript_host);

print 'Running on probe2trascript.pl on: '.`hostname`."\n";

my $transcript_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $transcript_host,
						       '-port'   => $transcript_port,
						       '-user'   => $transcript_user,
						       '-pass'   => $transcript_pass,
						       '-dbname' => $transcript_dbname);
my $oligo_db;

if ($oligo_host && $oligo_dbname && $oligo_user) {

  $oligo_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $oligo_host,
						 '-port'   => $oligo_port,
						 '-user'   => $oligo_user,
						 '-pass'   => $oligo_pass,
						 '-dbname' => $oligo_dbname);

} else {

  print "No oligo database specified, reading oligo features from $transcript_host:$transcript_port:$transcript_dbname\n";
  $oligo_db = $transcript_db;

}

my $xref_db;

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

check_names_match($oligo_db, $xref_db);

delete_existing_xrefs($oligo_db, $xref_db, @arrays) if ($delete);
delete_unmapped_entries($xref_db) if ($delete);

check_existing_and_exit($oligo_db, $xref_db, @arrays);

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

my %promiscuous_probesets;
my %transcripts_per_probeset;

my %transcript_ids;
my %transcript_probeset_count; # key: transcript:probeset value: count

my %arrays_per_probeset;
my %array_probeset_sizes;

my @unmapped_objects;

$| = 1; # auto flush stdout

my $i = 0;
my $last_pc = -1;

my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();

my @transcripts = @{$transcript_adaptor->fetch_all()};
my $total =  scalar(@transcripts);

open (LOG, ">${transcript_dbname}_probe2transcript.log");

print "Identified ".scalar(@transcripts)." for probe mappinng\n";
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

 
  foreach my $feature (@$oligo_features) {
    #next if ($transcript->strand() != $feature->strand()); # XXX log this

    my $probe = $feature->probe();
    my $probe_name = ${$probe->get_all_complete_names()}[0];
    my $probeset = $probe->probeset();

    my $transcript_key = $transcript->stable_id() . ":" . $probeset;

    my $probe_length = $probe->probelength();
    my $min_overlap = ($probe_length - $max_mismatches);

    my $exon_overlap = $rr->overlap_size("exonic", $feature->seq_region_start(), $feature->seq_region_end());
    my $utr_overlap  = $rr->overlap_size("utr",    $feature->seq_region_start(), $feature->seq_region_end());

    if ($exon_overlap >= $min_overlap) {

#	  print "Exon overlap $exon_overlap excedes minimum overlap $min_overlap\n";


      $transcript_probeset_count{$transcript_key}{$probe_name}++;

    } elsif ($utr_overlap >= $min_overlap) {

#	  print "UTR overlap $utr_overlap excedes minimum overlap $min_overlap\n";

      $transcript_probeset_count{$transcript_key}{$probe_name}++;

    } else { # must be intronic


#	  print "Intronic\n";

      print LOG "Unmapped intronic " . $transcript->stable_id . "\t" . $probeset . " probe length $probe_length\n";
      if (!$no_triage) {


		my $um_obj = new Bio::EnsEMBL::UnmappedObject(-type       => 'probe2transcript',
													  -analysis   => $analysis,
													  -identifier => $probeset,
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

# cache which arrays a probeset belongs to, and the sizes of probesets in different arrays
cache_arrays_per_probeset($oligo_db);


print "Writing xrefs\n";
my $um_cnt = 0;

# now loop over all the mappings and add xrefs for those that have a suitable number of matches
foreach my $key (keys %transcript_probeset_count) {

  my ($transcript, $probeset) = split (/:/, $key);

  # store one xref/object_xref for each array-probeset-transcript combination
  foreach my $array (split (/ /, $arrays_per_probeset{$probeset})) {

    next if (@arrays && find_in_list($array, @arrays,) == -1 );

    my $size_key = $array . ":" . $probeset;
    my $probeset_size = $array_probeset_sizes{$size_key};

    my $hits = scalar(keys %{$transcript_probeset_count{$key}});

    if ($hits / $probeset_size >= $mapping_threshold) {

      # only create xrefs for non-promiscuous probesets

      #XXX
      add_xref($transcript_ids{$transcript}, $probeset, $db_entry_adaptor, $array, $probeset_size, $hits);
      print LOG "$probeset\t$transcript\tmapped\t$probeset_size\t$hits\n";

     # if ($transcripts_per_probeset{$probeset} <= $max_transcripts_per_probeset) {
     #
     #   add_xref($transcript_ids{$transcript}, $probeset, $db_entry_adaptor, $array, $probeset_size, $hits);
     #   $transcripts_per_probeset{$probeset}++;
     #   print LOG "$probeset\t$transcript\tmapped\t$probeset_size\t$hits\n";
     #
     # } else {
     #
     #   print LOG "$probeset\t$transcript\tpromiscuous\t$probeset_size\t$hits\n";
     #   $promiscuous_probesets{$probeset} = $probeset;
     #   # TODO - remove mappings for probesets that end up being promiscuous
     #
     # }

      # TODO - write insufficient/promiscuous/orphan to unmapped_object ?

    } else {

      print LOG "$probeset\t$transcript\tinsufficient\t$probeset_size\t$hits\n";
       if (!$no_triage) {

		 my $um_obj = new Bio::EnsEMBL::UnmappedObject(-type       => 'probe2transcript',
													   -analysis   => $analysis,
													   -identifier => $probeset,
													   -summary    => "Insufficient hits",
													   -full_desc  => "Probe had an insufficient number of hits (probeset size = $probeset_size, hits = $hits)",
													   -ensembl_object_type => 'Transcript',
													   -ensembl_id => $transcript_ids{$transcript});
		 
		 #push @unmapped_objects, new Bio::EnsEMBL::UnmappedObject(-type       => 'probe2transcript',
		#						 -analysis   => $analysis,
		#						 -identifier => $probeset,
		#						 -summary    => "Insufficient hits",
		#						 -full_desc  => "Probe had an insufficient number of hits (probeset size = $probeset_size, hits = $hits)",
		#						 -ensembl_object_type => 'Transcript',
		#						 -ensembl_id => $transcript_ids{$transcript});


		 &cache_and_load_unmapped_objects($um_obj);

      }

    }

  }

}



#can we load first batch of unmapped_object here to save memmory

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
  my ($um_obj) = @_;

  push @unmapped_objects, $um_obj;

  if(scalar(@unmapped_objects) >10000){
	#print "Uploading " . scalar(@unmapped_objects) . " unmapped reasons to xref database\n";

	if($first_cache){
	  $unmapped_objects[0]->dbID('2000');
	  $first_cache = 0;
	}

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
		my $um_obj = new Bio::EnsEMBL::UnmappedObject(-type       => 'probe2transcript',
																 -analysis   => $analysis,
																 -identifier => $probeset,
																 -summary    => "No transcript mappings",
																 -full_desc  => "Probeset did not map to any transcripts");

		&cache_and_load_unmapped_objects($um_obj);

		#push @unmapped_objects, new Bio::EnsEMBL::UnmappedObject(-type       => 'probe2transcript',
		#														 -analysis   => $analysis,
		#														 -identifier => $probeset,
		#														 -summary    => "No transcript mappings",
		#														 -full_desc  => "Probeset did not map to any transcripts");
      }
    }
  }
}

# ----------------------------------------------------------------------

sub cache_arrays_per_probeset {

  my ($db) = @_;

  print "Caching arrays per probeset\n";

  #my $sth = $db->dbc()->prepare("SELECT op.probeset, oa.name, oa.probe_setsize FROM oligo_probe op, oligo_array oa WHERE oa.oligo_array_id=op.oligo_array_id GROUP BY op.probeset, oa.name ");

  #do not need distinct on count as we're linking by array?
  my $sth = $db->dbc()->prepare("SELECT op.probeset, oa.name, count(op.oligo_probe_id) FROM oligo_probe op, oligo_array oa WHERE oa.oligo_array_id=op.oligo_array_id GROUP BY op.probeset, oa.name");


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

# Delete existing xrefs & object xrefs. Use user-specified arrays if
# defined, otherwise all arrays.
# Assumes external_db.dbname == oligo_array.name

sub delete_existing_xrefs {

  my ($oligo_adaptor, $dbentry_adaptor, @arrays) = @_;

  # get array names if necessary
  my @arrays_to_delete;
  if (@arrays) {

    @arrays_to_delete = @arrays;

  } else {

    my $sth = $oligo_adaptor->dbc()->prepare("SELECT DISTINCT(name) FROM oligo_array ORDER BY name");
    $sth->execute();

    while(my @row = $sth->fetchrow_array()){

      push @arrays_to_delete, $row[0];

    }

    $sth->finish();

  }

  my $del_sth = $dbentry_adaptor->dbc()->prepare("DELETE x, ox FROM xref x, object_xref ox, external_db e WHERE x.xref_id=ox.xref_id AND e.external_db_id=x.external_db_id AND e.db_name = ?");

  foreach my $array (@arrays_to_delete) {

    print "Deleting xrefs and object_xrefs for $array\n";
    $del_sth->execute($array);

  }

  $del_sth->finish();

}

# ----------------------------------------------------------------------

# Delete existing entries in unmapped_object, unmapped_reason and analysis

sub delete_unmapped_entries {

  my ($xref_adaptor) = @_;


  #This deletes all unmapped objects, even if we're only (re-)running one array!

  my $del_sth = $xref_adaptor->dbc()->prepare("DELETE a, ur, uo FROM analysis a, unmapped_reason ur, unmapped_object uo WHERE a.logic_name = 'probe2transcript' AND a.analysis_id=uo.analysis_id AND uo.unmapped_reason_id=ur.unmapped_reason_id");

  print "Deleting unmapped records\n";
  $del_sth->execute();

  $del_sth->finish();

}

# ----------------------------------------------------------------------

# Check if there are already xrefs defined, and exit if there are.
# Use user-specified arrays if defined, otherwise all arrays.
# Assumes external_db.dbname == oligo_array.name

sub check_existing_and_exit {

  my ($oligo_adaptor, $dbentry_adaptor) = @_;

  # get array names if necessary
  my @arrays_to_check;
  if (@arrays) {

    @arrays_to_check = @arrays;

  } else {

    my $sth = $oligo_adaptor->dbc()->prepare("SELECT DISTINCT(name), probe_setsize FROM oligo_array");
    $sth->execute();

    while(my @row = $sth->fetchrow_array()){

      #this is not essential now as we're counting the size of each probeset
      warn "Array $row[0] does not have a probeset_size set, please rectify\n"  if($row[1] == 0);

      push @arrays_to_check, $row[0];

    }

    $sth->finish();

  }

  my $xref_sth = $dbentry_adaptor->dbc()->prepare("SELECT COUNT(*) FROM xref x, object_xref ox, external_db e WHERE x.xref_id=ox.xref_id AND e.external_db_id=x.external_db_id AND e.db_name = ?");

  foreach my $array (@arrays_to_check) {

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

# Check that the array names in oligo_db match the names in the external_db table in xref_db

sub check_names_match {

 my ($oligo_adaptor, $dbentry_adaptor) = @_;

 my $mismatch;

 # cache oligo_array.name entries
 my %oligo_array_names;

 my $oligo_sth = $oligo_adaptor->dbc()->prepare("SELECT DISTINCT(name) FROM oligo_array");
 $oligo_sth->execute();

 while(my @row = $oligo_sth->fetchrow_array()){

   $oligo_array_names{$row[0]} = $row[0];

 }

 $oligo_sth->finish();

 # and external_db names
 my %external_db_names;
 my $xref_sth = $dbentry_adaptor->dbc()->prepare("SELECT DISTINCT(db_name) FROM external_db");

 $xref_sth->execute();
 while (my @row2 = $xref_sth->fetchrow_array()) {

   $external_db_names{$row2[0]} = $row2[0];

 }

  $xref_sth->finish();

 # now compare them
 foreach my $oligo_array_name (keys %oligo_array_names) {
   if (!$external_db_names{$oligo_array_name}) {
     print "$oligo_array_name appears in oligo_array but not in external_db\n";
     $mismatch = 1;
   }
 }

 #foreach my $external_db_name (keys %external_db_names) {
 #  if (!$oligo_array_names{$external_db_name}) {
 #    print "$external_db_name appears in external_db but not in oligo_array\n";
 #  }
 #}

 if ($mismatch) {
   print "At least one oligo array name was found that does not appear in external_db; this needs to be rectified before this script can be run successfully.\n";
   exit(1);
 }

}

# ----------------------------------------------------------------------

# Remove mappings for a particular probeset

#sub remove_probeset_transcript_mappings {
#
#  my ($dba, $probeset) = @_;
#
#  my $p = 0;
#
#  my $sth = $dba->dbc()->prepare("DELETE FROM object_xref WHERE xref_id=? AND ensembl_object_type='Transcript' AND ensembl_id=?");
#
#  my $values = @{$dbentries_per_probeset{$probeset}};
#
#  foreach my $value (@{$dbentries_per_probeset{$probeset}}) {
#
#    my ($dbe_id, $transcript_id) = split(/:/, $value);
#
#    $sth->execute($dbe_id, $transcript_id);
#    $p++;
#
#  }
#
#  $sth->finish();
#
#}

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

  [--max_transcripts] Only use this many transcripts. Useful for debugging.

  [--no_triage]       Don't write to the unmapped_object/unmapped_reason tables.

  [--health_check]    Only do sanity checks, then stop.

  [--help]            This text.


EOF

  exit(0);

}

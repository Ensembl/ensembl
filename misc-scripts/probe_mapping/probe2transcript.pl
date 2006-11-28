use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

my ($transcript_host, $transcript_port, $transcript_user, $transcript_pass, $transcript_dbname,
    $oligo_host, $oligo_port, $oligo_user, $oligo_pass, $oligo_dbname,
    $xref_host, $xref_port, $xref_user, $xref_pass, $xref_dbname,
    $print, $max_mismatches, $utr_length, $max_probesets_per_transcript, $max_transcripts, @arrays, $delete);

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
           'utr_length=i'           => \$utr_length,
	   'max_probesets=i'        => \$max_probesets_per_transcript,
	   'max_transcripts=i'      => \$max_transcripts,
	   'arrays=s'               => \@arrays,
	   'print'                  => \$print,
	   'delete'                 => \$delete,
           'help'                   => sub { usage(); exit(0); });

$transcript_port ||= 3306; $oligo_port ||= 3306; $xref_port ||= 3306;

$max_mismatches ||= 1;

$utr_length ||= 2000;

$max_probesets_per_transcript ||= 100;

@arrays = split(/,/,join(',',@arrays));

usage() if(!$transcript_user || !$transcript_dbname || !$transcript_host);

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

delete_existing_xrefs($oligo_db, $xref_db, @arrays) if ($delete);

check_existing_and_exit($oligo_db, $xref_db, @arrays);

my $transcript_adaptor = $transcript_db->get_TranscriptAdaptor();
my $slice_adaptor = $transcript_db->get_SliceAdaptor();
my $oligo_feature_adaptor = $oligo_db->get_OligoFeatureAdaptor();
my $db_entry_adaptor = $xref_db->get_DBEntryAdaptor();

my %count;

my %promiscuous_probesets;
my %dbentries_per_probeset;

$| = 1; # auto flush stdout

my $i = 0;
my $last_pc = -1;

my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();

my @transcripts = @{$transcript_adaptor->fetch_all()};
my $total =  scalar(@transcripts);

print "Mapping, percentage complete: ";

foreach my $transcript (@transcripts) {

  my $pc = int ((100 * $i) / $total);

  if ($pc > $last_pc) {
    print "$pc ";
    $last_pc = $pc;
  }

  last if ($max_transcripts && $i >= $max_transcripts);

  my $stable_id = $transcript->stable_id();

  my $slice = $transcript->feature_Slice(); # note not ->slice() as this gets whole chromosome!
  my $extended_slice = $slice->expand(0, $utr_length);

  my $exons = $transcript->get_all_Exons();

  $rr->flush();

  foreach my $exon (@$exons) {
    my $start = $exon->seq_region_start() - 13; # XXX ceil(25/2) like Craig's - do this properly
    my $end = $exon->seq_region_end() + 13;
    $rr->check_and_register("exonic", $start, $end);
  }

  $rr->check_and_register("utr", $extended_slice->end()-$utr_length, $extended_slice->end());

  my $oligo_features;
  if (@arrays) {
    $oligo_features = $extended_slice->get_all_OligoFeatures(@arrays);
  } else {
    $oligo_features = $extended_slice->get_all_OligoFeatures();
  }

  foreach my $feature (@$oligo_features) {

    my $probeset = $feature->probeset();

    next if ($promiscuous_probesets{$probeset});

    my $probe_length = $feature->probe()->probelength();
    my $min_overlap = ($probe_length - $max_mismatches);

    my $exon_overlap = $rr->overlap_size("exonic", $feature->seq_region_start(), $feature->seq_region_end());
    my $utr_overlap  = $rr->overlap_size("utr",    $feature->seq_region_start(), $feature->seq_region_end());

    if ($exon_overlap >= $min_overlap) {

      $count{'exonic'}++;
      add_xref($transcript, $feature, $db_entry_adaptor) if (!$print);

    } elsif ($utr_overlap >= $min_overlap) {

      $count{'utr'}++;
      add_xref($transcript, $feature, $db_entry_adaptor) if (!$print);

    } else { # must be intronic

      $count{'intronic'}++;

    }
  }

  # TODO - make external_db array names == array names in OligoArray!
  # change oligo_array.name to be oligo_array.external_db_id ?

  $i++;

}

print "\n";

print_stats();

# ----------------------------------------------------------------------

sub print_stats {

  my $e = $count{'exonic'};
  my $i = $count{'intronic'};
  my $u = $count{'utr'};
  my $t = $e + $i + $u;

  print "Total probesets: $t\n\n";
  print "Exonic:  \t$e\t" . pc($e, $t) . "%\n";
  print "Intronic:\t$i\t" . pc($i, $t) . "%\n";
  print "UTR:     \t$u\t" . pc($u, $t) . "%\n\n";

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

  my ($transcript, $feature, $dbea) = @_;

  my $probeset = $feature->probeset();

  # store one xref/object_xref for each array-probeset-transcript combination

  foreach my $array (@{$feature->probe()->get_all_Arrays()}) {

    my $array_name = $array->name();
    next if (@arrays && find_in_list($array_name, @arrays,) == -1 );

    # TODO - this only works if external_db db_name == array name; currently needs
    # some manual hacking to acheive this

    my $dbe = new Bio::EnsEMBL::DBEntry( -adaptor              => $dbea,
					 -primary_id           => $probeset,
					 -version              => "1",
					 -dbname               => $array_name,
					 -release              => "1",
					 -display_id           => $probeset,
					 -description          => undef,
					 -primary_id_linkable  => 1,
					 -display_id_linkable  => 0,
					 -priority             => 1,
					 -db_display_name      => $array_name,
					 -info_type            => "MISC",  # TODO - change to PROBE when available
					 -info_text            => "probe2transcript.pl test");


    $dbea->store($dbe, $transcript->dbID(), "Transcript");

    # store the dbID of the newly created DBEntry in %dbentries_per_probeset
    # so that promiscuous ones can be removed later; note format of value is
    # $dbe->dbID:$transcript->dbID
    push @{$dbentries_per_probeset{$probeset}}, $dbe->dbID() . ":" . $transcript->dbID();

    # if any probesets map to more than 100 transcripts, ignore them in future and
    # delete existing mappings to them
    if (scalar(@{$dbentries_per_probeset{$probeset}}) > $max_probesets_per_transcript) {
      $promiscuous_probesets{$probeset} = $probeset;
      #remove_probeset_transcript_mappings($t_db, $probeset);
    }
  }

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

    my $sth = $oligo_adaptor->dbc()->prepare("SELECT DISTINCT(name) FROM oligo_array");
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

    my $sth = $oligo_adaptor->dbc()->prepare("SELECT DISTINCT(name) FROM oligo_array");
    $sth->execute();

    while(my @row = $sth->fetchrow_array()){

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

# Remove mappings for a particular probeset

sub remove_probeset_transcript_mappings {

  my ($dba, $probeset) = @_;

  my $p = 0;

  my $sth = $dba->dbc()->prepare("DELETE FROM object_xref WHERE xref_id=? AND ensembl_object_type='Transcript' AND ensembl_id=?");

  my $values = @{$dbentries_per_probeset{$probeset}};

  foreach my $value (@{$dbentries_per_probeset{$probeset}}) {

    my ($dbe_id, $transcript_id) = split(/:/, $value);

    $sth->execute($dbe_id, $transcript_id);
    $p++;

  }

  $sth->finish();

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

  GENERAL MAPPING OPTIONS:

  [--mismatches]      Allow up to this number of mismatches, inclusive.
                      Defaults to 1.
		
  [--utr_length]      Search this many bases downstream of the transcript
                      coding region as well. Defaults to 2000.
		
  [--max_probesets]   Don't store mappings to any 'promiscuous' probesets that map
                      to more than this number of transcripts. Defaults to 100.

  [--arrays]          Comma-separated list of arrays to use. Defaults to all arrays.

  MISCELLANEOUS:

  [--delete]          Delete existing xrefs and object_xrefs. No deletion is done by default.

  [--print]           Print information about mapping, don't store in database.

  [--max_transcripts] Only use this many transcripts. Useful for debugging.

  [--help]            This text.


EOF

  exit(0);

}

use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

my ($s_host, $s_port, $s_user, $s_pass, $s_dbname, $t_host, $t_port, $t_user, $t_pass, $t_dbname, $print, $max_mismatches, $utr_length, $max_probesets_per_transcript, $max_transcripts, @arrays);

GetOptions('source_host=s'      => \$s_host,
           'source_user=s'      => \$s_user,
           'source_port=i'      => \$s_port,
           'source_pass=s'      => \$s_pass,
           'source_dbname=s'    => \$s_dbname,
	   'target_host=s'      => \$t_host,
           'target_user=s'      => \$t_user,
           'target_port=i'      => \$t_port,
           'target_pass=s'      => \$t_pass,
           'target_dbname=s'    => \$t_dbname,
	   'mismatches=i'       => \$max_mismatches,
           'utr_length=i'       => \$utr_length,
	   'max_probesets=i'    => \$max_probesets_per_transcript,
	   'max_transcripts=i'  => \$max_transcripts,
	   'arrays=s'           => \@arrays,
	   'print'              => \$print,
           'help'               => sub { usage(); exit(0); });

$s_port ||= 3306; $t_port ||= 3306;

$max_mismatches ||= 1;

$utr_length ||= 2000;

$max_probesets_per_transcript ||= 100;

@arrays = split(/,/,join(',',@arrays));

usage() if(!$s_user || !$s_dbname || !$s_host);


my $s_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $s_host,
					      '-port'   => $s_port,
					      '-user'   => $s_user,
					      '-pass'   => $s_pass,
					      '-dbname' => $s_dbname);

my $t_db;

if ($t_host && $t_dbname && $t_user) {

  $t_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $t_host,
					     '-port'   => $t_port,
					     '-user'   => $t_user,
					     '-pass'   => $t_pass,
					     '-dbname' => $t_dbname);

} else {

  print "No target database specified, writing xrefs to $s_host:$s_port:$s_dbname\n";
  $t_db = $s_db;

}

my $transcript_adaptor = $t_db->get_TranscriptAdaptor();
my $slice_adaptor = $t_db->get_SliceAdaptor();
my $oligo_feature_adaptor = $s_db->get_OligoFeatureAdaptor();
my $db_entry_adaptor = $t_db->get_DBEntryAdaptor();

my %exonic_probesets;
my %utr_probesets;
my %intronic_probesets;

my %promiscuous_probesets;
my %dbentries_per_probeset;

my $i = 0;

my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();

foreach my $transcript (@{$transcript_adaptor->fetch_all()}) {

  print "$i\n" if ($i % 1000 == 0);

  last if ($max_transcripts && $i >= $max_transcripts);

  my $stable_id = $transcript->stable_id();

  my $slice = $transcript->feature_Slice(); # note not ->slice() as this gets whole chromosome!
  my $extended_slice = $slice->expand(0, $utr_length);

  my $exons = $transcript->get_all_Exons();

  $rr->flush();

  foreach my $exon (@$exons) {
    $rr->check_and_register("exonic", $exon->seq_region_start, $exon->seq_region_end);
  }

  $rr->check_and_register("utr", $extended_slice->end()-$utr_length, $extended_slice->end());

  my $oligo_features = $extended_slice->get_all_OligoFeatures();

  foreach my $feature (@$oligo_features) {

    my $probeset = $feature->probeset();

    next if ($promiscuous_probesets{$probeset});

    my $probe_length = $feature->probe()->probelength();
    my $min_overlap = ($probe_length - $max_mismatches);

    my $exon_overlap = $rr->overlap_size("exonic", $feature->seq_region_start(), $feature->seq_region_end());
    my $utr_overlap = $rr->overlap_size("utr", $feature->seq_region_start(), $feature->seq_region_end());

    if ($exon_overlap >= $min_overlap) {

      $exonic_probesets{$probeset} = $stable_id;
      add_xref($transcript, $feature, $db_entry_adaptor) if (!$print);

    } elsif ($utr_overlap > $min_overlap) {

      $utr_probesets{$probeset} = $stable_id;
      add_xref($transcript, $feature, $db_entry_adaptor) if (!$print);

    } else { # must be intronic

      $intronic_probesets{$probeset} = $stable_id;

    }
  }

  # TODO - make external_db array names == array names in OligoArray!
  # change oligo_array.name to be oligo_array.external_db_id ?

  $i++;

}

print_stats();

# ----------------------------------------------------------------------

sub print_stats {

  my $e = scalar(keys(%exonic_probesets));
  my $i = scalar(keys(%intronic_probesets));
  my $u = scalar(keys(%utr_probesets));
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

    next if (@arrays && find_in_list($array->name(), @arrays,) == -1 );

    # TODO - get db name from array name; for now just use AFFY_HG_U133A

    my $dbe = new Bio::EnsEMBL::DBEntry( -adaptor              => $dbea,
					 -primary_id           => $probeset,
					 -version              => "1",
					 -dbname               => "AFFY_HG_U133A", # XXXXX use proper array name
					 -release              => "1",
					 -display_id           => $probeset,
					 -description          => undef,
					 -primary_id_linkable  => 1,
					 -display_id_linkable  => 0,
					 -priority             => 1,
					 -db_display_name      => "AFFY_HG_U133A", # XXXXX
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
      remove_probeset_transcript_mappings($t_db, $probeset);
    }
  }

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

  READING OLIGO FEATURES:

   --source_host          The database server to read oligo features from.
		
   [--source_port]        The port to use for reading features. Defaults to 3306.
		
   --source_user          Database username for reading features.
		
   --source_pass          Password for source_user, if required.
		
   --source_dbname        Database name to read features from.

  WRITING XREFS:

   --target_host          The database server to write xrefs to.
		
   [--target_port]        The port to use for writing xrefs.. Defaults to 3306.
		
   --target_user          Database username for xrefs. Must allow writing.
		
   --target_pass          Password for target_user, if required.
		
   --target_dbname        Database name to write xrefs to.

  Note that if no target_host etc is specified, xrefs will be written to the
  database specified by the source_* parameters.

  GENERAL MAPPING OPTIONS:

  [--mismatches]      Allow up to this number of mismatches, inclusive.
                      Defaults to 1.
		
  [--utr_length]      Search this many bases downstream of the transcript
                      coding region as well. Defaults to 2000.
		
  [--max_probesets]   Don't store mappings to any 'promiscuous' probesets that map
                      to more than this number of transcripts. Defaults to 100.

  [--arrays]          Comma-separated list of arrays to use. Defaults to all arrays.

  MISCELLANEOUS:

  [--print]           Print information about mapping, don't store in database.

  [--max_transcripts] Only use this many transcripts. Useful for debugging.

  [--help]            This text.


EOF

  exit(0);

}

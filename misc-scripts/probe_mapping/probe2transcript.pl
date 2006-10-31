use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

my ($s_host, $s_port, $s_user, $s_pass, $s_dbname, $t_host, $t_port, $t_user, $t_pass, $t_dbname, $print, $max_mismatches, $utr_length);

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
	   'print'              => \$print,
           'help'               => sub { usage(); exit(0); });

$s_port ||= 3306; $t_port ||= 3306;

$max_mismatches ||= 1;

$utr_length ||= 2000;

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

my $i = 0;

foreach my $transcript (@{$transcript_adaptor->fetch_all()}) {

  print "$i\n" if ($i % 1000 == 0);

  my $stable_id = $transcript->stable_id();

  my $slice = $transcript->feature_Slice(); # note not ->slice() as this gets whole chromosome!
  my $extended_slice = $slice->expand(0, $utr_length);

  my $exons = $transcript->get_all_Exons();

  my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();

  foreach my $exon (@$exons) {
    $rr->check_and_register("exonic", $exon->seq_region_start, $exon->seq_region_end);
  }

  $rr->check_and_register("utr", $extended_slice->end()-$utr_length, $extended_slice->end());

  my $oligo_features = $extended_slice->get_all_OligoFeatures();

  foreach my $feature (@$oligo_features) {

    my $probeset = $feature->probeset();
    my $probe = $feature->probe();
    my $probe_length = $probe->probelength();
    my $min_overlap = ($probe_length - $max_mismatches);

    my $arrays = "";
    foreach my $array (@{$probe->get_all_Arrays()}) {
      $arrays .= $array->name() . " ";
    }

    my $exon_overlap = $rr->overlap_size("exonic", $feature->seq_region_start(), $feature->seq_region_end());
    my $utr_overlap = $rr->overlap_size("utr", $feature->seq_region_start(), $feature->seq_region_end());

    if ($exon_overlap >= $min_overlap) {

      $exonic_probesets{$probeset} = $stable_id;
      print join("\t", $stable_id, $probeset, $exon_overlap, $arrays, "exonic", "\n") if ($print);
      add_xref($transcript, $probe, $db_entry_adaptor) if (!$print);

    } elsif ($utr_overlap > $min_overlap) {

      $utr_probesets{$probeset} = $stable_id;
      print "$stable_id\t$probeset\t$utr_overlap\tutr\n" if ($print);
      add_xref($transcript, $probe, $db_entry_adaptor) if (!$print);

    } else { # must be intronic

      $intronic_probesets{$probeset} = $stable_id;
      print "$stable_id\t$probeset\t\tintronic\n" if ($print);

    }
  }

  # TODO - make external_db array names == array names in OligoArray!
  # change oligo_array.name to be oligo_array.external_db_id ?

  # TODO - remove 'promiscuous' probes

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

  my ($transcript, $probe, $dbea) = @_;

  # TODO - get db name from probe name; for now just use AFFY_HG_U133A

  my $dbe = new DBEntry( -adaptor              => $dbea,
			 -primary_id           => $probe->name(),
			 -version              => "1",
			 -dbname               => "AFFY_HG_U133A", # XXXXX
			 -release              => "1",
			 -display_id           => $probe->name(),
			 -description          => undef,
			 -primary_id_linkable  => 1,
			 -display_id_linkable  => 0,
			 -priority             => 1,
			 -db_display_name      => "AFFY_HG_U133A", # XXXXX
			 -info_type            => "MISC",  # TODO - change to PROBE when available
			 -info_text            => "");

  $dbea->store($dbe, $transcript->dbID(), "Transcript");

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

  Note that if no target_host etc is specified, xrefs will be written to the database specified
  by the source_* parameters.

  GENERAL OPTIONS:

  [--mismatches]   Allow up to this number of mismatches, inclusive. Defaults to 1.

  [--utr_length]   Search this many bases downstream of the transcript coding region as well. Defaults to 2000.

  [--print]        Print information about mapping, don't store in database.

  [--help]         This text.


EOF

  exit(0);

}

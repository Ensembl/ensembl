# Convert release 16-era schemas to use new non-clone/contig schema

use strict;
use warnings;

use DBI;
use Getopt::Long;

my ($host, $port, $user, $password, $source, $target, $verbose, $create, $clean, $check);
$host = "127.0.0.1";
$port = 5000;
$password = "";
$user = "ensro";

GetOptions ('host=s'      => \$host,
            'user=s'      => \$user,
            'password=s'  => \$password,
            'port=s'      => \$port,
            'source=s'    => \$source,
            'target=s'    => \$target,
            'verbose'     => \$verbose,
	    'clean'       => \$clean,
	    'check'       => \$check,
	    'create=s'    => \$create,
            'help'        => sub { &show_help(); exit 1;} );

die "Host must be specified"           unless $host;
die "Target schema must be specified"  unless $target;
die "Source schema be specifed"        unless $source;

# clean and create need to be done in a specific order
# they create their own db connections as necessary
&clean()  if $clean;
&create() if $create;

my $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$target", "$user", "$password") || die "Can't connect to target DB";
my $sth;

# ----------------------------------------------------------------------
# The coord_system table needs to be filled first

debug("Building coord_system table");
my @inserts = ('INSERT INTO coord_system (name, version, attrib) VALUES ("chromosome",  "NCBI33", "default_version,top_level")',
	       'INSERT INTO coord_system (name, version, attrib) VALUES ("supercontig", NULL,     "default_version")',
	       'INSERT INTO coord_system (name, version, attrib) VALUES ("clone",       NULL,     "default_version")',
	       'INSERT INTO coord_system (name, version, attrib) VALUES ("contig",      NULL,     "default_version,sequence")' );
foreach my $insert (@inserts) {
  execute($dbi, $insert);
}

# cache coord-system names to save lots of joins
debug("Caching coord_system IDs");
my %coord_system_ids; # hash with keys = coord-system names, values = internal IDs
$sth = $dbi->prepare("SELECT coord_system_id, name FROM $target.coord_system");
$sth->execute or die "Error when caching coord-system IDs";
while(my $row = $sth->fetchrow_hashref()) {
  my $id = $row->{"coord_system_id"};
  my $name = $row->{"name"};
  $coord_system_ids{$name} = $id;
}

# ----------------------------------------------------------------------
# Add info to the meta table that describes which features are in which
# co-ordinate systems

copy_table($dbi, "meta");
debug("Adding co-ordinate system information to meta table");
$sth = $dbi->prepare("INSERT INTO meta (meta_key, meta_value) VALUES ('feature.coord_system', ?)");
my @cs = ("gene|chromosome:NCBI33", "transcript|chromosome:NCBI33", "exon|chromosome:NCBI33", "dna_align_feature|contig", "marker_feature|contig", "simple_feature|contig", "protein_align_feature|contig", "repeat_feature|contig", "qtl_feature|chromosome");
foreach my $val (@cs) {
  $sth->execute($val);
}

# ----------------------------------------------------------------------
# Build the seq_region table
# This table is fundamental since many other tables reference it by seq_region_id

# For now we can just copy the contents of the contig table, with the appropriate
# coordinate system ID. Note contig IDs are copied as-is (i.e. NOT using autonumber)

debug("Copying source contig table to seq_region");
my $cs_id = $coord_system_ids{"contig"};
execute($dbi, "INSERT INTO seq_region SELECT contig_id, name, $cs_id, length from $source.contig");

# Similarly for the clone table - can use autonumber for the IDs as they're not referenced anywhere
$cs_id = $coord_system_ids{"clone"};
execute($dbi, "INSERT INTO seq_region (name, coord_system_id, length) SELECT CONCAT(cl.name, '.', cl.version), $cs_id, MAX(ctg.embl_offset)+ctg.length-1 FROM $source.clone cl, $source.contig ctg WHERE cl.clone_id=ctg.clone_id GROUP BY ctg.clone_id");

# And chromosomes
# Note old/new ID mapping is stored in %chromosme_id_old_new 
my %chromosome_id_old_new;
$cs_id = $coord_system_ids{"chromosome"};
$sth = $dbi->prepare("SELECT chromosome_id, name, length FROM $source.chromosome");
$sth->execute or die "Error when caching coord-system IDs";
while(my $row = $sth->fetchrow_hashref()) {
  my $old_id = $row->{"chromosome_id"};
  my $name   = $row->{"name"};
  my $length = $row->{"length"};
  execute($dbi, "INSERT INTO seq_region (name, coord_system_id, length) VALUES ('$name', $cs_id, $length)");
  my $new_id = $dbi->{'mysql_insertid'};
  $chromosome_id_old_new{$old_id} = $new_id;
}
#foreach my $id (keys (%chromosome_id_old_new)) {
#  print "OLD: $id   NEW: " . $chromosome_id_old_new{$id} . "\n";
#}

# ----------------------------------------------------------------------
# Gene
# Need to calculate start, end etc

debug("Building gene table");

my $sql =
  "INSERT INTO $target.gene " .
  "SELECT g.gene_id, g.type, g.analysis_id, e.contig_id, " .
  "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
  "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
  "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
  "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
  "       a.contig_ori*e.contig_strand as strand, " .
  "       g.display_xref_id " .
  "FROM   $source.transcript t, $source.exon_transcript et, $source.exon e, $source.assembly a, $source.gene g " .
  "WHERE  t.transcript_id = et.transcript_id " .
  "AND    et.exon_id = e.exon_id " .
  "AND    e.contig_id = a.contig_id " .
  "AND    g.gene_id = t.gene_id " . 
  "GROUP BY g.gene_id";
execute($dbi, $sql);

# ----------------------------------------------------------------------
# Transcript
# Need to calculate start, end etc

debug("Building transcript table");

$sql =
  "INSERT INTO $target.transcript " .
  "SELECT t.transcript_id, g.gene_id, e.contig_id, " .
  "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
  "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
  "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
  "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
  "       a.contig_ori*e.contig_strand as strand, " .
  "       g.display_xref_id " .
  "FROM   $source.transcript t, $source.exon_transcript et, $source.exon e, $source.assembly a, $source.gene g " .
  "WHERE  t.transcript_id = et.transcript_id " .
  "AND    et.exon_id = e.exon_id " .
  "AND    e.contig_id = a.contig_id " .
  "AND    g.gene_id = t.gene_id " .
  "GROUP BY t.transcript_id";
#print $sql . "\n";
execute($dbi, $sql);

# ----------------------------------------------------------------------
# Exon
# Translation to chromosomal co-ordinates should take care of sticky exons
# so new exon table will have <= number of rows as old exon table
debug("Building exon table");

$sql =
  "INSERT INTO $target.exon " .
  "SELECT e.exon_id, e.contig_id, " .
  "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
  "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
  "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
  "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
  "       a.contig_ori*e.contig_strand as strand, " .
  "       e.phase, e.end_phase " .
  "FROM   $source.transcript t, $source.exon_transcript et, $source.exon e, $source.assembly a, $source.gene g " .
  "WHERE  t.transcript_id = et.transcript_id " .
  "AND    et.exon_id = e.exon_id " .
  "AND    e.contig_id = a.contig_id " .
  "AND    g.gene_id = t.gene_id " .
  "GROUP BY e.exon_id";
#print $sql . "\n";
execute($dbi, $sql);

# ----------------------------------------------------------------------
# Translation
# Now includes transcript_id

debug("Building translation table");

$sql = 
  "INSERT INTO $target.translation " .
  "SELECT tl.translation_id, ts.transcript_id, tl.seq_start, tl.start_exon_id, tl.seq_end, tl.end_exon_id " .
  "FROM $source.transcript ts, $source.translation tl " .
  "WHERE ts.translation_id = tl.translation_id";
#print $sql . "\n";
execute($dbi, $sql);

# ----------------------------------------------------------------------
# Assembly

debug("Building assembly table");

# Need to do this row-by-row to get chromosome ID mapping
$sth = $dbi->prepare("SELECT * FROM $source.assembly");
$sth->execute or die "Error when reading assembly info from source DB";
while(my $row = $sth->fetchrow_hashref()) {

  $sql = "INSERT INTO $target.assembly VALUES (" . 
         $chromosome_id_old_new{$row->{"chromosome_id"}} . ", " .
	 $row->{"contig_id"} . ", " .
	 $row->{"chr_start"} . ", " .
	 $row->{"chr_end"} . ", " .
	 $row->{"contig_start"} . ", " .
	 $row->{"contig_end"} . ", " .
	 $row->{"contig_ori"} . ")";
  execute($dbi, $sql);

}

# ----------------------------------------------------------------------
# dna table

debug("Copying dna table");
execute($dbi, "INSERT INTO $target.dna SELECT dna_id, sequence FROM $source.dna");

# ----------------------------------------------------------------------
# Feature tables
# Note that we can just rename contig_* to set_region_* since the
# contig IDs were copied verbatim into seq_region

# simple_feature
debug("Copying simple_feature");
execute($dbi, "INSERT INTO $target.simple_feature (simple_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, display_label, analysis_id, score) SELECT simple_feature_id, contig_id, contig_start, contig_end, contig_strand, 'display_label', analysis_id, score FROM $source.simple_feature");

# repeat_feature
debug("Copying repeat_feature");
execute($dbi, "INSERT INTO $target.repeat_feature (repeat_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id, repeat_start, repeat_end, repeat_consensus_id, score) SELECT repeat_feature_id, contig_id, contig_start, contig_end, contig_strand, analysis_id, repeat_start, repeat_end, repeat_consensus_id, score FROM $source.repeat_feature");

# protein_align_feature
debug("Copying protein_align_feature");
execute($dbi, "INSERT INTO $target.protein_align_feature (protein_align_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, hit_name, cigar_line, evalue, perc_ident, score) SELECT protein_align_feature_id, contig_id, contig_start, contig_end, contig_strand, analysis_id, hit_start, hit_end, hit_name, cigar_line, evalue, perc_ident, score FROM $source.protein_align_feature");

# dna_align_feature
debug("Copying dna_align_feature");
execute($dbi, "INSERT INTO $target.dna_align_feature (dna_align_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, hit_name, hit_strand, cigar_line, evalue, perc_ident, score) SELECT dna_align_feature_id, contig_id, contig_start, contig_end, contig_strand, analysis_id, hit_start, hit_end, hit_name, hit_strand, cigar_line, evalue, perc_ident, score FROM $source.dna_align_feature");

# marker_feature
debug("Copying marker_feature");
execute($dbi, "INSERT INTO $target.marker_feature (marker_feature_id, marker_id, seq_region_id, seq_region_start, seq_region_end, analysis_id, map_weight) SELECT marker_feature_id, marker_id, contig_id, contig_start, contig_end, analysis_id, map_weight FROM $source.marker_feature");

# qtl_feature
# Note this uses the perviously constructed %chromosome_id_old_new hash for mapping
debug("Copying qtl_feature");
$sth = $dbi->prepare("SELECT * FROM $source.qtl_feature");
$sth->execute or die "Error when reading qtl_feature";
while(my $row = $sth->fetchrow_hashref()) {

  my $seq_region_id = $chromosome_id_old_new{$row->{"chromosome_id"}};
  execute($dbi, "INSERT INTO $target.qtl_feature(seq_region_id, start, end, qtl_id, analysis_id) VALUES (" . $seq_region_id . ", " . $row->{"start"} . ", " . $row->{"end"} . ", " . $row->{"qtl_id"} . ", " . $row->{"analysis_id"} . ")");

}

# ----------------------------------------------------------------------
# These tables now have seq_region_* instead of chromosome_*

debug("Translating karyotype");
$sth = $dbi->prepare("SELECT * FROM $source.karyotype");
$sth->execute or die "Error when reading karyotype";
while(my $row = $sth->fetchrow_hashref()) {
  my $seq_region_id = $chromosome_id_old_new{$row->{"chromosome_id"}};
  execute($dbi, "INSERT INTO $target.karyotype VALUES (" . $seq_region_id . ", " .
	         $row->{"chr_start"} . ", " .
	         $row->{"chr_end"} . ", " .
	         "'" . $row->{"band"} . "', " .
	         "'" . $row->{"stain"} . "')");
}

debug("Translating marker_map_location");
$sth = $dbi->prepare("SELECT * FROM $source.marker_map_location");
$sth->execute or die "Error when reading marker_map_location";
while(my $row = $sth->fetchrow_hashref()) {
  my $seq_region_id = $chromosome_id_old_new{$row->{"chromosome_id"}};
  execute($dbi, "INSERT INTO $target.marker_map_location VALUES (" .
	         $row->{"marker_id"} . ", " .
	         $row->{"map_id"} . ", " .
	         $seq_region_id . ", " .
	         $row->{"marker_synonym_id"} . ", " .
	         "'" . $row->{"position"} . "', " .
	         $row->{"lod_score"} . ")");
}

debug("Translating map_density");
$sth = $dbi->prepare("SELECT * FROM $source.map_density");
$sth->execute or die "Error when reading map_density";
while(my $row = $sth->fetchrow_hashref()) {
  my $seq_region_id = $chromosome_id_old_new{$row->{"chromosome_id"}};
  execute($dbi, "INSERT INTO $target.map_density VALUES (" .
	         $seq_region_id . ", " .
	         $row->{"chr_start"} . ", " .
	         $row->{"chr_end"} . ", " .
	         "'" . $row->{"type"} . "', " .
	         $row->{"value"} . ")");
}

# ----------------------------------------------------------------------
# These tables are copied as-is]
copy_table($dbi, "supporting_feature");
copy_table($dbi, "map");
copy_table($dbi, "analysis");
copy_table($dbi, "dnafrag");
copy_table($dbi, "exon_stable_id");
copy_table($dbi, "exon_transcript");
copy_table($dbi, "external_db");
copy_table($dbi, "external_synonym");
copy_table($dbi, "gene_archive");
copy_table($dbi, "gene_description");
copy_table($dbi, "gene_stable_id");
copy_table($dbi, "go_xref");
copy_table($dbi, "identity_xref");
copy_table($dbi, "interpro");
copy_table($dbi, "map");
copy_table($dbi, "mapannotation");
copy_table($dbi, "mapannotationtype");
copy_table($dbi, "mapfrag");
copy_table($dbi, "mapfrag_mapset");
copy_table($dbi, "mapping_session");
copy_table($dbi, "mapset");
copy_table($dbi, "marker");
copy_table($dbi, "marker_feature");
copy_table($dbi, "marker_synonym");
copy_table($dbi, "object_xref");
copy_table($dbi, "peptide_archive");
copy_table($dbi, "protein_feature");
copy_table($dbi, "qtl");
copy_table($dbi, "qtl_synonym");
copy_table($dbi, "repeat_consensus");
copy_table($dbi, "stable_id_event");
copy_table($dbi, "transcript_stable_id");
copy_table($dbi, "translation_stable_id");
copy_table($dbi, "xref");

# TODO finish

# ----------------------------------------------------------------------

&check() if $check;

$dbi->disconnect();

debug("Done");

# ----------------------------------------------------------------------
# Misc / utility functions

sub show_help {

  print "Usage: perl convert_seqstore.pl {options}\n";
  print "Where options are:\n";
  print "  --host {hostname} The database host.\n";
  print "  --user {username} The database user. Must have read permissions on both schema, and write permissions on the target schema\n";
  print "  --password {pass} The password for user, if required.\n";
  print "  --port {folder}   The database port to use.\n";
  print "  --source {schema} The name of the source schema\n";
  print "  --target {schema} The name of the target schema\n";
  print "  --clean           Remove target schema, which must have been specified with --target\n";
  print "  --create {file}   Create target schema, which must have been specified with --target, from SQL file\n";
  print "  --check           Check target schema for empty tables at end od run\n";
  print "  --verbose         Print extra output information\n";

}

# ----------------------------------------------------------------------

sub debug {

  my $str = shift;

  print $str . "\n" if $verbose;

}

# ----------------------------------------

sub execute {

  my ($dbcon, $sql) = @_;

  my $stmt = $dbcon->prepare($sql);
  $stmt->execute() || die "Unable to execute $sql";
  $stmt->finish();

}

# ----------------------------------------

sub copy_table {

  my ($dbcon, $table) = @_;

  debug("Copying $table with no modifications");
  my $sql = "INSERT INTO $target.$table SELECT * FROM $source.$table";
  execute($dbcon, $sql);

}

# ----------------------------------------

sub clean {

  if (!defined $target) {

    print "Target database must be specified with --target for --clean\n";
    exit(1);

  } else {

    debug("Removing $target");

    my $dbic = DBI->connect("dbi:mysql:host=$host;port=$port;", "$user", "$password") || die "Can't connect to DB";
    execute($dbic, "DROP DATABASE $target") || die "Error removing $target";
    $dbic->disconnect();

  }

}

# ----------------------------------------

sub create {

  if (!defined $target) {

    print "Target database must be specified with --target for --clean\n";
    exit(1);

  } else {

    debug("Creating database $target");

    my $dbic = DBI->connect("dbi:mysql:host=$host;port=$port;", "$user", "$password") || die "Can't connect to DB";
    execute($dbic, "CREATE DATABASE $target") || die "Error creating $target";
    $dbic->disconnect();

    debug("Building schema for $target from $create");
    die "Can't open $create" if (! -e $create);
    my $cmd = "mysql -u $user -p$password -h $host -P $port $target < $create";
    #debug($cmd);
    system ($cmd);

  }

}

# ----------------------------------------

sub check {

  my ($sth_source, $sth_target);
  my $dbi_source = DBI->connect("dbi:mysql:host=$host;port=$port;database=$source", "$user", "$password") || die "Can't connect to DB";
  my $dbi_target = DBI->connect("dbi:mysql:host=$host;port=$port;database=$target", "$user", "$password") || die "Can't connect to DB";

  $sth = $dbi_source->prepare("SHOW TABLES");
  $sth->execute or die "Error when listing tables";
  while(my @row = $sth->fetchrow_array()) {

    my $table_name = $row[0];

    $sth_source = $dbi_source->prepare("SELECT COUNT(*) FROM " . $table_name);
    $sth_source->execute or die "Error when counting rows in source " . $table_name;
    my @count = $sth_source->fetchrow_array();
    my $source_count = $count[0];

    $sth_target = $dbi_target->prepare("SELECT COUNT(*) FROM " . $table_name);
    my $res = $sth_target->execute();
    if (defined $res) {
      @count = $sth_target->fetchrow_array();
      my $target_count = $count[0];
      if ($source_count > $target_count) {
	print "Warning: " . $table_name . " has " . $source_count . " rows in " . $source . " but " . $target_count . " rows in " . $target . "\n";
      }
    }
  }

  $sth_source->finish();
  $sth_target->finish();
  $dbi_source->disconnect();
  $dbi_target->disconnect();

}

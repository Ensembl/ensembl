# Convert release 18-era schemas to use new non-clone/contig schema

use strict;
use warnings;
use DBI;

package SeqStoreConverter::BasicConverter;


sub new {
  my ( $class, $user, $pass, $host, $source, $target, $schema, $force, 
       $verbose ) = @_;

  my $self = bless {} $class;

  my $dbh = DBI->connect( "DBI:mysql:host=$host", $user, $password );

  $self->verbose( $verbose );
  $self->dbh( $dbh );
  $self->force( $force );
  $self->source( $source );
  $self->target( $target );
  $self->schema( $schema );
  $self->host( $host );
  $self->password( $pass);
  $self->user($user);


  #check to see if the destination and source databases exist already.
  my %dbs = map {$_->[0] => 1} @{$dbh->selectall_arrayref('show databases')};
  if( ! $dbs{$source} ) {
    die ("source db $source doesnt exist" );
  }

  if( $dbs{$target} ) {
    if( !$force ) {
      die( "target db $target already exists" );
    } else {
      $dbh->execute( "drop database $target" );
    }
  }
  
  create_target_database();
  return $self;
}


sub force {
  my $self = shift;
  $self->{'force'} = shift if(@_);
  return $self->{'force'};
}


sub dbh {
  my $self = shift;
  $self->{'dbh'} = shift if(@_);
  return $self->{'dbh'};
}

sub user {
  my $self = shift;
  $self->{'user'} = shift if(@_);
  return $self->{'user'};
}

sub host {
  my $self = shift;
  $self->{'host'} = shift if(@_);
  return $self->{'host'};
}

sub password {
  my $self = shift;
  $self->{'password'} = shift if(@_);
  return $self->{'password'};
}


sub verbose {
  my $self = shift;
  $self->{'verbose'} = shift if (@_);
  return $self->{'verbose'};
}

sub schema {
  my $self = shift;
  $self->{'schema'} = shift if (@_);
  return $self->{'schema'};
}

sub debug {
  my $self = shift;
  my $str = shift;

  print STDERR $str . "\n" if $self->verbose();
}


my %coord_system_ids; # hash with keys = coord-system names, values = internal IDs

# ----------------------------------------------------------------------
# The coord_system and meta_coord tables need to be filled first but
# how this is done varies from species to species

# copy the old meta table info first
copy_table($dbi, "meta");

build_species_coord_tables($species, $dbi, %coord_system_ids);

# ----------------------------------------------------------------------
# Build the seq_region table and the assembly table
# This table is fundamental since many other tables reference it by seq_region_id
# Exactly what is done is species-dependent

build_seq_region_and_assembly($species, $dbi, %coord_system_ids);

#-----------------------------------------------------------------------------------------------
# seq_region_attrib
# Place the clones HTG phase into seq_region attrib.
# Presently nothing else needs to go in here

# only done for certain species as this only makes sense if there are "real" clones
# done by checking if there is a "clone" coordinate system in the /target/ database


debug( "Translating mappannotationtype" );
#
# first copy the old map annotation types into the attrib type table
#
execute( $dbi,
         "INSERT INTO $target.attrib_type( attrib_type_id, code, name, description ) " .
         "SELECT mapannotationtype_id, code, name, description " .
         "FROM $source.mapannotationtype " );

debug( "Populating seq_region_attrib table" );

if (species_has_clone_table($species)) {
  debug("Species has clone information, transferring HTG phase data");
  
  my $sth = $dbi->prepare("SELECT distinct(htg_phase) FROM clone");
  $sth->execute();

  # now add a new attrib type for HTG phase
  execute( $dbi,
           "INSERT INTO $target.attrib_type( code, name, description ) " .
           " VALUES ('htg_phase', 'HTG Phase', 'High Throughput Genome Phase')");
  
  #insert the htg phase values for the clones into seq_region_attrib
  execute( $dbi,
           "INSERT INTO $target.seq_region_attrib( seq_region_id, attrib_type_id, value) " .
           "SELECT tmp_cln.new_id, attrib_type.attrib_type_id, cln.htg_phase " .
           "FROM   $target.tmp_cln_map tmp_cln, $target.attrib_type attrib_type, $source.clone cln " .
           "WHERE  cln.clone_id = tmp_cln.old_id " .
           "AND    attrib_type.code = 'htg_phase'");
} else {
  debug("Species has no clone information, not transferring HTG phase data");
}

# ----------------------------------------------------------------------
# Gene
# Need to calculate start, end etc

debug("Building gene table");

my $sql =
  "INSERT INTO $target.gene " .
  "SELECT g.gene_id, g.type, g.analysis_id, tcm.new_id, " .
  "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
  "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
  "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
  "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
  "       a.contig_ori*e.contig_strand as strand, " .
  "       g.display_xref_id " .
  "FROM   $source.transcript t, $source.exon_transcript et, $source.exon e, $source.assembly a, $source.gene g, $target.tmp_chr_map tcm " .
  "WHERE  t.transcript_id = et.transcript_id " .
  "AND    et.exon_id = e.exon_id " .
  "AND    e.contig_id = a.contig_id " .
  "AND    g.gene_id = t.gene_id ";

if ($species ne "fugu") {
  $sql .= "AND    a.chromosome_id = tcm.old_id ";
}

$sql .= "GROUP BY g.gene_id";

execute($dbi, $sql);

# ----------------------------------------------------------------------
# Transcript
# Need to calculate start, end etc

debug("Building transcript table");

$sql =
  "INSERT INTO $target.transcript " .
  "SELECT t.transcript_id, t.gene_id, tcm.new_id, " .
  "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
  "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
  "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
  "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
  "       a.contig_ori*e.contig_strand as strand, " .
  "       t.display_xref_id " .
  "FROM   $source.transcript t, $source.exon_transcript et, $source.exon e, $source.assembly a, $target.tmp_chr_map tcm " .
  "WHERE  t.transcript_id = et.transcript_id " .
  "AND    et.exon_id = e.exon_id " .
  "AND    e.contig_id = a.contig_id ";

if ($species ne "fugu") {
  $sql .= "AND    a.chromosome_id = tcm.old_id ";
}

$sql .=  "GROUP BY t.transcript_id";

execute($dbi, $sql);

# ----------------------------------------------------------------------
# Exon
# Translation to chromosomal co-ordinates should take care of sticky exons
# so new exon table will have <= number of rows as old exon table
debug("Building exon table");

$sql =
  "INSERT INTO $target.exon " .
  "SELECT e.exon_id, tcm.new_id, " .
  "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
  "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
  "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
  "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
  "       a.contig_ori*e.contig_strand as strand, " .
  "       e.phase, e.end_phase " .
  "FROM   $source.transcript t, $source.exon_transcript et, $source.exon e, $source.assembly a, $source.gene g, $target.tmp_chr_map tcm " .
  "WHERE  t.transcript_id = et.transcript_id " .
  "AND    et.exon_id = e.exon_id " .
  "AND    e.contig_id = a.contig_id " .
  "AND    g.gene_id = t.gene_id ";

if ($species ne "fugu") {
  $sql .= "AND    a.chromosome_id = tcm.old_id ";
}
$sql .= "GROUP BY e.exon_id";

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
# dna table

debug("Translating dna");
execute($dbi, "INSERT INTO $target.dna " .
	"SELECT c.contig_id as seq_region_id, d.sequence as sequence " .
	"FROM   $source.dna d, $source.contig c " .
	"WHERE  c.dna_id = d.dna_id");

# ----------------------------------------------------------------------
# Feature tables
# Note that we can just rename contig_* to set_region_* since the
# contig IDs were copied verbatim into seq_region

# For some reason mysql refuses to use the index on large tables sometimes
# following copies like the following.
# So: drop the indexes first and then re-add them after
# It is probably faster this way anyway


# simple_feature
debug("Translating simple_feature");
execute($dbi, "INSERT INTO $target.simple_feature (simple_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, display_label, analysis_id, score) SELECT simple_feature_id, contig_id, contig_start, contig_end, contig_strand, display_label, analysis_id, score FROM $source.simple_feature");

# repeat_feature
if (!$noindex) {
  debug("Dropping indexes on repeat_feature");
  execute($dbi, "ALTER TABLE $target.repeat_feature DROP INDEX seq_region_idx");
  execute($dbi, "ALTER TABLE $target.repeat_feature DROP INDEX repeat_idx");
  execute($dbi, "ALTER TABLE $target.repeat_feature DROP INDEX analysis_idx");
}

debug("Translating repeat_feature");
execute($dbi, "INSERT INTO $target.repeat_feature (repeat_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id, repeat_start, repeat_end, repeat_consensus_id, score) SELECT repeat_feature_id, contig_id, contig_start, contig_end, contig_strand, analysis_id, repeat_start, repeat_end, repeat_consensus_id, score FROM $source.repeat_feature");

if (!$noindex) {
  debug("Readding indexes on repeat_feature");
  execute($dbi, "ALTER TABLE $target.repeat_feature ADD INDEX seq_region_idx( seq_region_id, seq_region_start) ");
  execute($dbi, "ALTER TABLE $target.repeat_feature ADD INDEX repeat_idx( repeat_consensus_id )");
  execute($dbi, "ALTER TABLE $target.repeat_feature ADD INDEX analysis_idx(analysis_id)");
}

# protein_align_feature
if (!$noindex) {
  debug("Dropping indexes on protein_align_feature");
  execute($dbi, "ALTER TABLE $target.protein_align_feature DROP INDEX seq_region_idx");
  execute($dbi, "ALTER TABLE $target.protein_align_feature DROP INDEX hit_idx");
  execute($dbi, "ALTER TABLE $target.protein_align_feature DROP INDEX ana_idx");
  execute($dbi, "ALTER TABLE $target.protein_align_feature DROP INDEX score_idx");
}

debug("Translating protein_align_feature");
execute($dbi, "INSERT INTO $target.protein_align_feature (protein_align_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, hit_name, cigar_line, evalue, perc_ident, score) SELECT protein_align_feature_id, contig_id, contig_start, contig_end, contig_strand, analysis_id, hit_start, hit_end, hit_name, cigar_line, evalue, perc_ident, score FROM $source.protein_align_feature");

if (!$noindex) {
  debug("Readding indexes on protein_align_feature");
  execute($dbi, "ALTER TABLE $target.protein_align_feature ADD index seq_region_idx(seq_region_id, seq_region_start)");
  execute($dbi, "ALTER TABLE $target.protein_align_feature ADD index hit_idx(hit_name)");
  execute($dbi, "ALTER TABLE $target.protein_align_feature ADD index ana_idx(analysis_id)");
  execute($dbi, "ALTER TABLE $target.protein_align_feature ADD index score_idx(score)");
}

# dna_align_feature

if (!$noindex) {
  debug("Dropping indexes on dna_align_feature");
  execute($dbi, "ALTER TABLE $target.dna_align_feature DROP INDEX seq_region_idx");
  execute($dbi, "ALTER TABLE $target.dna_align_feature DROP INDEX hit_idx");
  execute($dbi, "ALTER TABLE $target.dna_align_feature DROP INDEX ana_idx");
  execute($dbi, "ALTER TABLE $target.dna_align_feature DROP INDEX score_idx");
}

debug("Translating dna_align_feature");
execute($dbi, "INSERT INTO $target.dna_align_feature (dna_align_feature_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, hit_name, hit_strand, cigar_line, evalue, perc_ident, score) SELECT dna_align_feature_id, contig_id, contig_start, contig_end, contig_strand, analysis_id, hit_start, hit_end, hit_name, hit_strand, cigar_line, evalue, perc_ident, score FROM $source.dna_align_feature");

if (!$noindex) {
  debug("Readding indexes on dna_align_feature");
  execute($dbi, "ALTER TABLE $target.dna_align_feature ADD index seq_region_idx(seq_region_id, seq_region_start)");
  execute($dbi, "ALTER TABLE $target.dna_align_feature ADD index hit_idx(hit_name)");
  execute($dbi, "ALTER TABLE $target.dna_align_feature ADD index ana_idx(analysis_id)");
  execute($dbi, "ALTER TABLE $target.dna_align_feature ADD index score_idx(score)");
}

# marker_feature
debug("Translating marker_feature");
execute($dbi, "INSERT INTO $target.marker_feature (marker_feature_id, marker_id, seq_region_id, seq_region_start, seq_region_end, analysis_id, map_weight) SELECT marker_feature_id, marker_id, contig_id, contig_start, contig_end, analysis_id, map_weight FROM $source.marker_feature");

# qtl_feature
# Note this uses chromosome coords so we have to join with tmp_chr_map to get the mapping
debug("Translating qtl_feature");

execute($dbi,
	"INSERT INTO $target.qtl_feature(seq_region_id, seq_region_start, seq_region_end, qtl_id, analysis_id) " .
	"SELECT tcm.new_id, " .
	"       q.start, q.end, q.qtl_id, q.analysis_id " .
	"FROM $target.tmp_chr_map tcm, $source.qtl_feature q " .
	"WHERE tcm.old_id = q.chromosome_id");

# ----------------------------------------------------------------------
# These tables now have seq_region_* instead of chromosome_*

debug("Translating karyotype");
execute($dbi,
	"INSERT INTO $target.karyotype " .
	"SELECT null, tcm.new_id, " .
	"       k.chr_start, k.chr_end, k.band, k.stain " .
	"FROM $target.tmp_chr_map tcm, $source.karyotype k " .
	"WHERE tcm.old_id = k.chromosome_id");


debug("Translating marker_map_location");
execute($dbi,
	"INSERT INTO $target.marker_map_location " .
	"SELECT mml.marker_id, mml.map_id, " .
	"       c.name, " .
	"       mml.marker_synonym_id, mml.position, mml.lod_score " .
	"FROM $source.chromosome c, $source.marker_map_location mml " .
	"WHERE c.chromosome_id = mml.chromosome_id");

debug("Translating map_density");
execute($dbi,
	"INSERT INTO $target.map_density " .
	"SELECT tcm.new_id, ".
	"       md.chr_start, md.chr_end, md.type, md.value " .
	"FROM $target.tmp_chr_map tcm, $source.map_density md " .
	"WHERE tcm.old_id = md.chromosome_id");


debug( "Translating mapfrag" );
execute( $dbi,
	 "INSERT INTO $target.misc_feature( misc_feature_id, seq_region_id, " . 
	 "            seq_region_start, seq_region_end, seq_region_strand ) " .
	 "SELECT m.mapfrag_id, sr.seq_region_id, m.seq_start, m.seq_end, m.orientation " .
	 "FROM   $source.mapfrag m, $target.seq_region sr, $source.dnafrag d " . 
	 "WHERE  m.dnafrag_id = d.dnafrag_id " .
	 "AND    d.name = sr.name " );

debug( "Translating mapset" );
execute( $dbi,
	 "INSERT INTO $target.misc_set( misc_set_id, code, name, description, " .
	 "                              max_length ) " .
	 "SELECT mapset_id, code, name, description, max_length " . 
	 "FROM $source.mapset ms " );

debug( "Translating mapannotation" );
execute( $dbi,
	 "INSERT INTO $target.misc_attrib( misc_feature_id, attrib_type_id, " . 
	 "                    value ) ". 
	 "SELECT mapfrag_id, mapannotationtype_id, value " .
	 "FROM $source.mapannotation" );

debug( "Translating mapfrag_mapset" );
execute( $dbi,
	 "INSERT INTO $target.misc_feature_misc_set( misc_feature_id, misc_set_id ) ".
	 "SELECT mapfrag_id, mapset_id ".
	 "FROM $source.mapfrag_mapset " );

# ----------------------------------------------------------------------
# prediction_transcript / prediction_exon
# prediction transcripts for anopheles are in contig coords; need to convert to chromosomal
debug( "Translating prediction_exon" );

if ($species eq "anopheles") {

  # TODO

} else {

  execute( $dbi, "INSERT INTO prediction_exon ".
	   "( prediction_transcript_id, seq_region_id, seq_region_start, seq_region_end, " .
	   "  seq_region_strand, start_phase, score, p_value, exon_rank ) " .
	   "SELECT prediction_transcript_id, contig_id, contig_start, contig_end, contig_strand, " .
	   "       start_phase, score, p_value, exon_rank " .
	   "FROM $source.prediction_transcript" );

}

debug("Translating prediction_transcript");

if ($species eq "anopheles") {

  # TODO

} else {

  execute( $dbi, "INSERT INTO prediction_transcript ".
	   "( prediction_transcript_id, seq_region_id, seq_region_start, seq_region_end, " .
	   " seq_region_strand, analysis_id ) " .
	   " SELECT prediction_transcript_id, contig_id, MIN(contig_start), MAX(contig_end), contig_strand, analysis_id ".
	   "FROM $source.prediction_transcript " .
	   "GROUP BY prediction_transcript_id ");

}

#-----------------------------------------------------------------
# remove the unused created and modified dates from the stable ids

execute( $dbi, "INSERT INTO exon_stable_id " .
	 " (exon_id, stable_id, version) " .
	 "SELECT exon_id, stable_id, version " .
	 "FROM $source.exon_stable_id" );


execute( $dbi, "INSERT INTO gene_stable_id " .
	 " (gene_id, stable_id, version) " .
	 "SELECT gene_id, stable_id, version " .
	 "FROM $source.gene_stable_id" );


# ----------------------------------------------------------------------
# These tables are copied as-is


copy_table($dbi, "supporting_feature");
copy_table($dbi, "map");
copy_table($dbi, "analysis");
copy_table($dbi, "exon_transcript");
copy_table($dbi, "external_db");
copy_table($dbi, "external_synonym");
copy_table($dbi, "gene_archive");
copy_table($dbi, "gene_description");
copy_table($dbi, "go_xref");
copy_table($dbi, "identity_xref");
copy_table($dbi, "interpro");
copy_table($dbi, "map");
copy_table($dbi, "mapping_session");
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

# ----------------------------------------------------------------------

&check() if $check;

$dbi->disconnect();

debug("Done");

# ----------------------------------------------------------------------


# ----------------------------------------

sub execute {

  my ($dbcon, $sql) = @_;

  my $stmt;

  eval {
    $stmt = $dbcon->prepare($sql);
    $stmt->execute();
    $stmt->finish();
  };

  if ($@) {
    my ($p,$l,$f) = caller;
    warn("WARNING: Unable to execute [$sql]:\n" .
         "         $@" .
         "         Called from $f line $l\n" );
    return 0;
  }
  return 1;
}

# ----------------------------------------

sub copy_table {
  my ($dbcon, $table) = @_;

  debug("Copying $table");
  my $sql = "INSERT INTO $target.$table SELECT * FROM $source.$table";
  execute($dbcon, $sql) or warn("WARNING: Could not copy table $table\n");
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

sub create_target_database {
  my $self = shift;
  $self->dbh()->execute( "create database ".$self->target() );



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
  while (my @row = $sth->fetchrow_array()) {

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

# ----------------------------------------------------------------------
# Species-specific operations

# Build coord_system and meta_coord tables for different species

sub build_species_coord_tables {

  my ($species, $dbi) = @_;

  $species = lc($species);

  # get default assembly from meta table
  my $ass_def;
  my $stmt = $dbi->prepare("SELECT meta_value FROM $source.meta WHERE meta_key='assembly.default'");
  my $res = $stmt->execute();
  if (defined $res) {
    my @row = ($stmt->fetchrow_array());
    $ass_def = $row[0];
    debug("Assembly default for $species: $ass_def\n");
  }

  if (!defined $ass_def || $ass_def eq "") {
    warn("Cannot get assembly.default from meta table for $species");
  }

  my (@coords, %cs, @assembly_mappings);

  if ($species eq "human" ) {

    @coords = ('("chromosome","' . $ass_def . '","default_version,top_level")',
	       '("supercontig", NULL,          "default_version")',
	       '("clone",       NULL,          "default_version")',
	       '("contig",      NULL,          "default_version,sequence_level")' );

    @assembly_mappings =  ("chromosome:$ass_def|contig",
			   "clone|contig",
			   "supercontig|contig");

  } elsif ($species eq "rat") {

    @coords = ('("chromosome","' . $ass_def . '","default_version,top_level")',
	       '("supercontig", NULL,          "default_version")',
	       '("contig",      NULL,          "default_version,sequence_level")' );

    @assembly_mappings =  ("chromosome:$ass_def|contig",
			   "clone|contig",
			   "supercontig|contig");

  } elsif ($species eq "mouse") {

    @coords = ('("chromosome","' . $ass_def . '","default_version,top_level")',
	       '("supercontig", NULL,          "default_version")',
	       '("clone",       NULL,          "default_version")',
	       '("contig",      NULL,          "default_version,sequence_level")' );

    @assembly_mappings =  ("chromosome:$ass_def|contig",
			   "clone|contig",
			   "supercontig|contig");

  } elsif ($species eq "fugu") {

    @coords = ('("scaffold","' . $ass_def .   '","default_version,top_level")');

    @assembly_mappings =  ();

  } elsif ($species eq "anopheles") {

    @coords = ('("chromosome","' . $ass_def  .'","default_version,top_level")',
	       '("scaffold",    NULL,          "default_version")',
	       '("chunk",       NULL,          "default_version,sequence_level")');

    @assembly_mappings =  ("chromosome:$ass_def|scaffold",
			   "scaffold|chunk");

  } elsif ($species eq "zebrafish") {

    @coords = ('("scaffold","' . $ass_def .   '","default_version,top_level,sequence_level")');

    @assembly_mappings =  ();

  } elsif ($species eq "elegans") {

    @coords = ('("chromosome","' . $ass_def . '","default_version,top_level")',
	       '("contig",      NULL,         "default_version,sequence_level")' );

    @assembly_mappings =  ("chromosome:$ass_def|contig");

  } elsif ($species eq "briggsae") {

    @coords = ('("supercontig","' . $ass_def . '","default_version,top_level")',
	       '("contig",       NULL,          "default_version,sequence_level")' );

    @assembly_mappings =  ("supercontig|contig");

  } elsif ($species eq "drosophila") {

    @coords = ('("chromosome","' . $ass_def . '","default_version,top_level")',
	       '("chunk",       NULL,          "default_version,sequence_level")' );

    @assembly_mappings =  ("chromosome|chunk");


  } else {

    warn("\nWARNING: species-specific settings not yet defined for $species - can't build coord_system table!\n\n");

  }

  # ----------------------------------------
  # Meta_coord table

  if ($species eq "human" || $species eq "mouse" || $species eq "rat" || $species eq "drosophila") {

    %cs = (gene                  => 'chromosome',
	   transcript            => 'chromosome',
	   exon               	 => 'chromosome',
	   dna_align_feature     => 'contig',
	   protein_align_feature => 'contig',
	   marker_feature        => 'contig',
	   simple_feature        => 'contig',
	   repeat_feature        => 'contig',
	   qtl_feature           => 'chromosome',
	   misc_feature          => 'chromosome',
	   prediction_transcript => 'contig',
	   karyotype             => 'chromosome');

   } elsif ($species eq "fugu") {

     %cs = (gene                  => 'scaffold',
	    transcript            => 'scaffold',
	    exon               	 => 'scaffold',
	    dna_align_feature     => 'scaffold',
	    protein_align_feature => 'scaffold',
	    marker_feature        => 'scaffold',
	    simple_feature        => 'scaffold',
	    repeat_feature        => 'scaffold',
	    qtl_feature           => 'scaffold',
	    misc_feature          => 'scaffold',
	    prediction_transcript => 'scaffold',
	    karyotype             => 'scaffold');

   } elsif ($species eq "anopheles") {

     %cs = (gene                  => 'chromosome',
	    transcript            => 'chromosome',
	    exon                  => 'chromosome',
	    dna_align_feature     => 'chunk',
	    protein_align_feature => 'chunk',
	    marker_feature        => 'chunk',
	    simple_feature        => 'chunk',
	    repeat_feature        => 'chunk',
	    qtl_feature           => 'chromosome',
	    misc_feature          => 'chromosome',
	    prediction_transcript => 'chunk',
	    karyotype             => 'chromosome');

   } elsif ($species eq "zebrafish") {

     %cs = (gene                  => 'scaffold',
	    transcript            => 'scaffold',
	    exon               	 => 'scaffold',
	    dna_align_feature     => 'scaffold',
	    protein_align_feature => 'scaffold',
	    marker_feature        => 'scaffold',
	    simple_feature        => 'scaffold',
	    repeat_feature        => 'scaffold',
	    qtl_feature           => 'scaffold',
	    misc_feature          => 'scaffold',
	    prediction_transcript => 'scaffold',
	    karyotype             => 'scaffold');

   } elsif ($species eq "elegans") {

     %cs = (gene                  => 'chromosome',
	    transcript            => 'chromosome',
	    exon               	 => 'chromosome',
	    dna_align_feature     => 'contig',
	    protein_align_feature => 'contig',
	    marker_feature        => 'contig',
	    simple_feature        => 'contig',
	    repeat_feature        => 'contig',
	    qtl_feature           => 'chromosome',
	    misc_feature          => 'chromosome',
	    prediction_transcript => 'contig',
	    karyotype             => 'chromosome');

   } elsif ($species eq "briggsae") {

     %cs = (gene                  => 'supercontig',
	    transcript            => 'supercontig',
	    exon               	  => 'supercontig',
	    dna_align_feature     => 'contig',
	    protein_align_feature => 'contig',
	    marker_feature        => 'contig',
	    simple_feature        => 'contig',
	    repeat_feature        => 'contig',
	    qtl_feature           => 'supercontig',
	    misc_feature          => 'supercontig',
	    prediction_transcript => 'contig',
	    karyotype             => 'supercontig');

   } else {

    warn("\nWARNING: species-specific settings not yet defined for $species - can't build meta_coord table!\n\n");

  }

  # ----------------------------------------

  debug("Building coord_system table for " . $species);
  foreach my $coord (@coords) {
    execute($dbi, 'INSERT INTO coord_system (name, version, attrib) VALUES ' . $coord);
  }

  # ----------------------------------------------------------------------
  # cache coord-system names to save lots of joins
  debug("Caching coord_system IDs");
  $sth = $dbi->prepare("SELECT coord_system_id, name FROM $target.coord_system");
  $sth->execute or die "Error when caching coord-system IDs";
  while (my $row = $sth->fetchrow_hashref()) {
    my $id = $row->{"coord_system_id"};
    my $name = $row->{"name"};
    $coord_system_ids{$name} = $id;
  }

  debug("Building meta_coord table for " . $species);
  $sth = $dbi->prepare("INSERT INTO $target.meta_coord VALUES (?, ?)");
  foreach my $val (keys %cs) {
    $sth->execute($val, $coord_system_ids{$cs{$val}});
  }

  debug("Building adding assembly.mapping entries to meta table for " . $species);
  foreach my $mapping (@assembly_mappings) {
    execute($dbi, "INSERT INTO $target.meta(meta_key, meta_value) VALUES (\"assembly.mapping\", \"$mapping\")");
  }

}

# ----------------------------------------------------------------------

sub build_seq_region_and_assembly() {

  my ($species, $dbi) = @_;

  if ($species eq "human" ) {

    # Human is straight coord system -> coord system conversion
    # since system was originally designed for human assembly:
    # contig      -> contig
    # clone       -> clone
    # supercontig -> supercontig
    # chromosome  -> chromosome

    chromosome_to_seq_region();
    supercontig_to_seq_region();
    clone_to_seq_region();
    contig_to_seq_region();

    assembly_contig_chromosome();
    assembly_contig_clone();
    assembly_contig_supercontig();

  } elsif ($species eq "rat") {

    # TODO
    warn("seq_region and assembly tables not built for rat!");

  } elsif ($species eq "mouse") {

    # TODO
    warn("seq_region and assembly tables not built for mouse!");

  } elsif ($species eq "fugu") {

    # Fugu has only a single coordinate system - the scaffold
    # contig      -> scaffold
    # clone       -> x
    # supercontig -> x
    # chromosome  -> x

    contig_to_seq_region('scaffold');

    # fugu has no assemblies so don't need to do anything else

  } elsif ($species eq "anopheles") {

    # contig      -> chunk
    # clone       -> scaffold
    # supercontig -> x
    # chromosome  -> chromosome

    chromosome_to_seq_region();
    contig_to_seq_region("chunk");
    clone_to_seq_region("scaffold");

    assembly_contig_chromosome();
    assembly_contig_clone();

  } elsif ($species eq "zebrafish") {

    # contig      -> ???
    # clone       -> ???
    # supercontig -> ???
    # chromosome  -> ???

    # TODO

  } elsif ($species eq "elegans") {
    
    # contig      -> x
    # clone       -> clone
    # supercontig -> x
    # chromosome  -> chromosome
    
    chromosome_to_seq_region();
    clone_to_seq_region();

    assembly_contig_chromosome();

  } elsif ($species eq "briggsae") {

    supercontig_to_seq_region();
    contig_to_seq_region();

    assembly_contig_supercontig();

  } elsif ($species eq "drosophila") {

    chromosome_to_seq_region();
    contig_to_seq_region("chunk");

    assembly_contig_chromosome();

  } else {

    warn("\nWARNING: species-specific settings not yet defined for $species - can't build seq_region / assembly !\n\n");

  }

  # Note use the following to check if the assembly table looks sane:
  # SELECT a.asm_seq_region_id, a.cmp_seq_region_id, a.asm_start, a.asm_end, a.cmp_start, a.cmp_end, s.name, s.length, cs.name AS coord_system FROM glenn_new_schema.assembly a, glenn_new_schema.seq_region s, glenn_new_schema.coord_system cs WHERE a.asm_seq_region_id=s.seq_region_id AND s.coord_system_id=cs.coord_system_id;

}

# ----------------------------------------------------------------------

sub chromosome_to_seq_region() {

  debug("Transforming chromosome table into seq_region");

  # target coord_system will have a different ID
  my $target_cs_name = shift || "chromosome";
  my $cs_id = get_id_for_coord_system($target_cs_name);

  # Note old/new ID mapping is stored in %chromosme_id_old_new; it turns out it is vastly
  # quicker to store the mappings in a temporary table and join to it rather than
  # doing a row-by-row INSERT using this hash

  my %chromosome_id_old_new;
  die "Error getting coord_system_id for chromosome" if !defined $cs_id || $cs_id eq "";
  my $sth = $dbi->prepare("SELECT chromosome_id, name, length FROM $source.chromosome");
  $sth->execute or die "Error when building chromosome ID/ seq_region ID map";
  while (my $row = $sth->fetchrow_hashref()) {
    my $old_id = $row->{"chromosome_id"};
    my $name   = $row->{"name"};
    my $length = $row->{"length"};
    execute($dbi, "INSERT INTO seq_region (name, coord_system_id, length) VALUES ('$name', $cs_id, $length)");
    my $new_id = $dbi->{'mysql_insertid'};
    $chromosome_id_old_new{$old_id} = $new_id;
  }

  # store this hash in a temporary table to save having to do row-by-row inserts later
  my $create_sql = 
    "CREATE TEMPORARY TABLE $target.tmp_chr_map (" .
      "old_id INT, new_id INT,".
	"INDEX new_idx (new_id))";
  $sth = $dbi->prepare($create_sql);
  $sth->execute or die "Error when creating temporary chr_map table";
  $sth = $dbi->prepare("INSERT INTO $target.tmp_chr_map (old_id, new_id) VALUES (?, ?)");
  while (my ($old_id, $new_id) = each %chromosome_id_old_new) {
    $sth->execute($old_id, $new_id) || die "Error writing to tmp_chr_map";
  }

}

# ----------------------------------------------------------------------

sub supercontig_to_seq_region() {

  # target coord_system will have a different ID
  my $target_cs_name = shift || "supercontig";
  my $cs_id = get_id_for_coord_system($target_cs_name);

  # Supercontigs - need to store new (seq_region) ID->name mapping for later
  debug("Transforming supercontigs into seq_region");
  my %superctg_name_id;
  $sth = $dbi->prepare("SELECT superctg_name, " .
		       "MAX(superctg_end)-MIN(superctg_start)+1 AS length " .
		       "FROM $source.assembly " .
		       "GROUP BY superctg_name");
  $sth->execute or die "Error when building supercontig name / seq region ID map";
  while (my $row = $sth->fetchrow_hashref()) {
    my $name = $row->{"superctg_name"};
    my $length = $row->{"length"};
    execute($dbi, "INSERT INTO seq_region (name, coord_system_id, length) VALUES ('$name', $cs_id, $length)");
    my $new_id = $dbi->{'mysql_insertid'};
    $superctg_name_id{$name} = $new_id;
  }

  # store this hash in a temporary table to save having to do row-by-row inserts later
  my $create_sql =
    "CREATE TEMPORARY TABLE $target.tmp_superctg_map (" .
      "name VARCHAR(255), new_id INT, ".
	"INDEX new_idx (new_id))";
  $sth = $dbi->prepare($create_sql);
  $sth->execute or die "Error when creating temporary tmp_superctg_map table";
  $sth = $dbi->prepare("INSERT INTO $target.tmp_superctg_map (name, new_id) VALUES (?, ?)");
  while (my ($name, $new_id) = each %superctg_name_id) {
    $sth->execute($name, $new_id) || die "Error writing to tmp_superctg_map";
    #print "$name\t$new_id\n";
  }

}

# ----------------------------------------------------------------------

sub clone_to_seq_region() {

  debug("Transforming clone table into seq_region");

  # target coord_system will have a different ID
  my $target_cs_name = shift || "clone";
  my $cs_id = get_id_for_coord_system($target_cs_name);

  my %clone_id_old_new;
  die "Error getting coord_system_id for clone" if !defined $cs_id || $cs_id eq "";
  $sth = $dbi->prepare
    ("SELECT cl.clone_id,
             CONCAT(cl.embl_acc, '.', cl.embl_version),
             MAX(ctg.embl_offset)+ctg.length-1,
     FROM   $source.clone cl, $source.contig ctg
		 WHERE  cl.clone_id = ctg.clone_id GROUP BY ctg.clone_id");
  $sth->execute or die "Error when building clone ID / seq_region ID map";

  my ($clone_id, $embl_acc, $name, $length);
  $sth->bind_columns(\$clone_id, \$embl_acc, \$name, \$length);

  while ($sth->fetch()) {
    execute($dbi, "INSERT INTO seq_region (name, coord_system_id, length) VALUES ('$embl_acc', $cs_id, $length)");
    my $new_id = $dbi->{'mysql_insertid'};
    $clone_id_old_new{$clone_id} = $new_id;
  }

  # store this hash in a temporary table to save having to do row-by-row inserts later
  my $create_sql =
    "CREATE TEMPORARY TABLE $target.tmp_cln_map (" .
      "old_id INT, new_id INT, ".
	"INDEX new_idx (new_id))";
  $sth = $dbi->prepare($create_sql);
  $sth->execute or die "Error when creating temporary tmp_cln_map table";
  $sth = $dbi->prepare("INSERT INTO $target.tmp_cln_map (old_id, new_id) VALUES (?, ?)");
  while (my ($old_id, $new_id) = each %clone_id_old_new) {
    $sth->execute($old_id, $new_id) || die "Error writing to tmp_cln_map";
    #print "$old_id\t$new_id\n";
  }

}

# ----------------------------------------------------------------------

sub contig_to_seq_region() {

  debug("Copying source contig table to seq_region");

  # target coord_system will have a different ID
  my $target_cs_name = shift || "contig";

  my $cs_id = get_id_for_coord_system($target_cs_name);

  die "Error getting coord_system_id for contig" if !defined $cs_id || $cs_id eq "";
  execute($dbi, "INSERT INTO seq_region SELECT contig_id, name, $cs_id, length from $source.contig");

}

# ----------------------------------------------------------------------

sub assembly_contig_chromosome() {

  debug("Building assembly table - contig/chromosome");

  execute($dbi,
	  "INSERT INTO $target.assembly " .
	  "SELECT tcm.new_id, " .	# asm_seq_region_id (old-new chromosome ID mapping)
	  "a.contig_id, " .	# cmp_seq_region_id
	  "a.chr_start, " .	# asm_start
	  "a.chr_end, " .		# asm_end
	  "a.contig_start, " .	# cmp_start
	  "a.contig_end, " .	# cmp_end
	  "a.contig_ori " .	# ori
	  "FROM $target.tmp_chr_map tcm, $source.assembly a, $source.contig c " .
	  "WHERE tcm.old_id = a.chromosome_id " .
	  "AND c.contig_id = a.contig_id "); # only copy assembly entries that refer to valid contigs

}

# ----------------------------------------------------------------------

sub assembly_contig_clone() {

  debug("Building assembly table - contig/clone");

  execute($dbi,
	  "INSERT INTO $target.assembly " .
	  "SELECT tcm.new_id, " .	# asm_seq_region_id (old-new clone ID mapping)
	  "ctg.contig_id, ".	# cmp_seq_region_id
	  "ctg.embl_offset, " .	# asm_start
	  "ctg.embl_offset+ctg.length-1, " . # asm_end
	  "1, " .			# cmp_start
	  "ctg.length, " .	# cmp_end
	  "1 " .	# ori - contig always positively oriented on the clone
	  "FROM $target.tmp_cln_map tcm, " .
	  "$source.clone cln, $source.contig ctg " .
	  "WHERE tcm.old_id = cln.clone_id " .
	  "AND cln.clone_id = ctg.clone_id");

}

# ----------------------------------------------------------------------

sub assembly_contig_supercontig() {

  debug("Building assembly table - contig/supercontig");

  execute($dbi,
	  "INSERT INTO $target.assembly " .
	  "SELECT tsm.new_id, " .	# asm_seq_region_id (superctg name-seq_region_id mapping)
	  "a.contig_id, " .	# cmp_seq_region_id
	  "a.superctg_start, " .	# asm_start
	  "a.superctg_end, " .	# asm_end
	  "a.contig_start, " .	# cmp_start
	  "a.contig_end, " .	# cmp_end
	  "a.contig_ori " .	# ori
	  "FROM $target.tmp_superctg_map tsm, $source.assembly a, $source.contig c " .
	  "WHERE tsm.name = a.superctg_name " .
	  "AND c.contig_id = a.contig_id "); # only copy assembly entries that refer to valid contigs

}

# ----------------------------------------------------------------------

sub rename_coord_system() {

  my ($old, $new) = @_;

  execute($dbi, "UPDATE coord_system SET name='" . $new . "' WHERE name='" . $old . "'");

}

# ----------------------------------------------------------------------

sub get_id_for_coord_system() {

  my $name = shift;

  my $sth = $dbi->prepare("SELECT coord_system_id FROM coord_system WHERE name='" . $name . "'");
  $sth->execute or die "Cannot get coord_system_id for $name";
  my @row = $sth->fetchrow_array();
  return $row[0];

}

# ----------------------------------------------------------------------

sub species_has_clone_table() {

  my $species = shift;

  my $sth = $dbi->prepare("SELECT * FROM coord_system WHERE name='clone'");
  $sth->execute or die "Cannot check for existence of clone coordinate system in target schema";
  if ($sth->fetchrow_array()) {
    return 1;
  }
  return 0;

}

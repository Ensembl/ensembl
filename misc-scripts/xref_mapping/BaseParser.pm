package BaseParser;

use strict;

use DBI;
use Digest::MD5 qw(md5_hex);
use File::Path;
use POSIX qw(strftime);

use SwissProtParser;
use RefSeqParser;
use RefSeqGPFFParser;

my $host = "ecs1g";
my $port = 3306;
my $database = "glenn_test_xref";
my $user = "ensadmin";
my $password = "ensembl";

my $base_dir = ".";

my $dbi;
my %dependent_sources;
my %taxonomy2species_id;
my %name2species_id;

run() if (!defined(caller()));

# --------------------------------------------------------------------------------
# Get info about files to be parsed from the database

sub run {

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT s.source_id, su.source_url_id, s.name, su.url, su.checksum, su.parser FROM source s, source_url su WHERE s.download='Y' AND su.source_id=s.source_id ORDER BY s.name");
  $sth->execute();
  my ($source_id, $source_url_id, $name, $url, $checksum, $parser);
  $sth->bind_columns(\$source_id, \$source_url_id, \$name, \$url, \$checksum, \$parser);
  my $last_type = "";
  my $dir;
  while (my @row = $sth->fetchrow_array()) {

    # Download each source into the appropriate directory for parsing later
    # Also delete previous working directory if we're starting a new source type
    my $type = $name;

    $dir = $base_dir . "/" . $type;
    rmtree $dir if ($type ne $last_type);
    $last_type = $type;
    mkdir $dir if (!-e $dir);

    my ($file) = $url =~ /.*\/(.*)/;

    print "Downloading $url to $dir/$file\n";
    my $result = system("wget", "--quiet","--directory-prefix=$dir", $url);

    # if the file is compressed, the FTP server may or may not have automatically uncompressed it
    # TODO - read .gz file directly? open (FILE, "zcat $file|") or Compress::Zlib
    if ($file =~ /(.*)\.gz$/) {
      print "Uncompressing $dir/$file\n";
      system("gunzip", "$dir/$file");
      $file = $1;
    }

    # compare checksums and parse/upload if necessary
    # need to check file size as some .SPC files can be of zero length
    my $file_cs = md5sum("$dir/$file");
    if (!defined $checksum || $checksum ne $file_cs) {

      if (-s "$dir/$file") {

	print "Checksum for $file does not match, parsing\n";

	update_source($dbi, $source_url_id, $file_cs, $file);

	print "Parsing $file with $parser\n";
	$parser->run("$dir/$file", $source_id);

      } else {
	
	print $file . " has zero length, skipping\n";

      }

    } else {

      print $file . " has not changed, skipping\n";

    }
  }

  # remove last working directory
  # TODO reinstate after debugging
  #rmtree $dir;

}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "BaseParser";

  return $self;

}

# --------------------------------------------------------------------------------
# Get source ID for a particular file; matches url field 

sub get_source_id_for_filename {

  my ($self, $file) = @_;

  my $sql = "SELECT s.source_id FROM source s, source_url su WHERE su.source_id=s.source_id AND su.url LIKE  '%/" . $file . "%'"; 
  #print $sql . "\n";
  my $sth = dbi()->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $source_id;
  if (defined @row) {
    $source_id = $row[0];
  } else {
    warn("Couldn't get source ID for file $file\n");
    $source_id = -1;
  }

  return $source_id;

}

# --------------------------------------------------------------------------------
# Upload xrefs to the database

sub upload_xrefs {

  my ($self, @xrefs) = @_;

  my $dbi = dbi();

  if ($#xrefs > -1) {
    # remove all existing xrefs with same source ID

    # TODO re-instate deletion


    my $source_id = $xrefs[0]->{SOURCE_ID};
    my $del_sth = $dbi->prepare("DELETE FROM xref WHERE source_id=$source_id");
    #$del_sth->execute();

    # upload new ones
    my $xref_sth = $dbi->prepare("INSERT INTO xref (accession,label,description,source_id,species_id) VALUES(?,?,?,?,?)");
    my $pri_insert_sth = $dbi->prepare("INSERT INTO primary_xref VALUES(?,?,?,?,?)");
    my $pri_update_sth = $dbi->prepare("UPDATE primary_xref SET sequence=? WHERE xref_id=?");
    my $syn_sth = $dbi->prepare("INSERT INTO synonym VALUES(?,?,?)");
    my $dep_sth = $dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");

    local $xref_sth->{RaiseError}; # disable error handling here as we'll do it ourselves
    local $xref_sth->{PrintError};

    foreach my $xref (@xrefs) {

      # Create entry in xref table and note ID
      $xref_sth->execute($xref->{ACCESSION},
			 $xref->{LABEL},
			 $xref->{DESCRIPTION},
			 $xref->{SOURCE_ID},
			 $xref->{SPECIES_ID});

      # If there was an error, an xref with the same acc & source already exists.
      # If so, find its ID, otherwise get ID of xref just inserted
      my $xref_id = insert_or_select($xref_sth, $dbi->err, $xref->{ACCESSION}, $xref->{SOURCE_ID});

      # create entry in primary_xref table with sequence; if this is a "cumulative"
      # entry it may already exist, and require an UPDATE rather than an INSERT
      if (primary_xref_id_exists($xref_id)) {
	
	$pri_update_sth->execute($xref->{SEQUENCE}, $xref_id) || die $dbi->errstr;
	
      } else {
	
	$pri_insert_sth->execute($xref_id,
				 $xref->{SEQUENCE},
				 $xref->{SEQUENCE_TYPE},
				 $xref->{STATUS},
				 $xref->{SOURCE_ID}) || die $dbi->errstr;
      }

      # if there are synonyms, create xrefs for them and entries in the synonym table
      foreach my $syn (@{$xref->{SYNONYMS}}) {

	$xref_sth->execute($syn,
			   "",
			   "",
			   $xref->{SOURCE_ID},
			   $xref->{SPECIES_ID});

	my $syn_xref_id = insert_or_select($xref_sth, $dbi->err, $syn, $xref->{SOURCE_ID});

	$syn_sth->execute($xref_id, $syn_xref_id, $xref->{SOURCE_ID} ) || die $dbi->errstr;

      }				# foreach syn

      # if there are dependent xrefs, add xrefs and dependent xrefs for them
      foreach my $depref (@{$xref->{DEPENDENT_XREFS}}) {

	my %dep = %$depref;

	$xref_sth->execute($dep{ACCESSION},
			   "",
			   "",
			   $dep{SOURCE_ID},
			   $xref->{SPECIES_ID});

	my $dep_xref_id = insert_or_select($xref_sth, $dbi->err, $dep{ACCESSION}, $dep{SOURCE_ID});
			
	$dep_sth->execute($xref_id, $dep_xref_id, '', $dep{SOURCE_ID} ) || die $dbi->errstr;
	# TODO linkage anntation?

      }				# foreach dep


    }				# foreach xref

    $del_sth->finish() if defined $del_sth;
    $xref_sth->finish() if defined $xref_sth;
    $pri_insert_sth->finish() if defined $pri_insert_sth;
    $pri_update_sth->finish() if defined $pri_update_sth;

  }

}

# --------------------------------------------------------------------------------
# Get & cache a hash of all the source names for dependent xrefs (those that are
# in the source table but don't have an associated URL etc)

sub get_dependent_xref_sources {

  my $self = shift;

  if (!defined %dependent_sources) {

    my $dbi = dbi();
    my $sth = $dbi->prepare("SELECT name,source_id FROM source WHERE download='N'");
    $sth->execute() || die $dbi->errstr;
    while(my @row = $sth->fetchrow_array()) {
      my $source_name = $row[0];
      my $source_id = $row[1];
      $dependent_sources{$source_name} = $source_id;
    }
  }

  return %dependent_sources;

}

# --------------------------------------------------------------------------------
# Get & cache a hash of all the species IDs & taxonomy IDs.

sub taxonomy2species_id {

  my $self = shift;

  if (!defined %taxonomy2species_id) {

    my $dbi = dbi();
    my $sth = $dbi->prepare("SELECT species_id, taxonomy_id FROM species");
    $sth->execute() || die $dbi->errstr;
    while(my @row = $sth->fetchrow_array()) {
      my $species_id = $row[0];
      my $taxonomy_id = $row[1];
      $taxonomy2species_id{$taxonomy_id} = $species_id;
    }
  }

  return %taxonomy2species_id;

}

# --------------------------------------------------------------------------------
# Get & cache a hash of all the species IDs & species names.

sub name2species_id {

  my $self = shift;

  if (!defined %name2species_id) {

    my $dbi = dbi();
    my $sth = $dbi->prepare("SELECT species_id, name FROM species");
    $sth->execute() || die $dbi->errstr;
    while(my @row = $sth->fetchrow_array()) {
      my $species_id = $row[0];
      my $name = $row[1];
      $name2species_id{$name} = $species_id;
    }
  }

  return %name2species_id;

}

# --------------------------------------------------------------------------------
# Update a row in the source table

sub update_source {

  my ($dbi, $source_url_id, $checksum, $file) = @_;
  open(FILE, $file);
  my $file_date = POSIX::strftime('%Y%m%d%H%M%S', localtime((stat($file))[9]));
  close(FILE);

  my $sql = "UPDATE source_url SET checksum='" . $checksum . "', file_modified_date='" . $file_date . "', upload_date=NOW() WHERE source_url_id=" . $source_url_id;
  # TODO release?

  $dbi->prepare($sql)->execute() || die $dbi->errstr;

}


# --------------------------------------------------------------------------------

sub dbi {

  my $self = shift;

  if (!defined $dbi) {
    $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$database",
			"$user",
			"$password",
			{'RaiseError' => 1}) || die "Can't connect to database";
  }

  return $dbi;

}

# --------------------------------------------------------------------------------

sub md5sum {

  my $file = shift;
  open(FILE, $file);
  binmode(FILE);
  my $md5 = Digest::MD5->new->addfile(*FILE)->hexdigest();
  close(FILE);

  return $md5;

}

# --------------------------------------------------------------------------------

sub get_xref_id_by_accession_and_source {

  my ($acc, $source_id) = @_;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=?");
  $sth->execute($acc, $source_id) || die $dbi->errstr;
  my @row = $sth->fetchrow_array();
  my $xref_id = $row[0];

  return $xref_id;

}

# --------------------------------------------------------------------------------
# If there was an error, an xref with the same acc & source already exists.
# If so, find its ID, otherwise get ID of xref just inserted

sub insert_or_select {

  my ($sth, $error, $acc, $source) = @_;

  my $id;

  if ($error) {

    $id = get_xref_id_by_accession_and_source($acc, $source);
    #print "Got existing xref id " . $id . " for " . $acc . " " . $source . "\n";

  } else {
	
    $id = $sth->{'mysql_insertid'};

  }

  return $id;

}

# --------------------------------------------------------------------------------

sub primary_xref_id_exists {

  my $xref_id = shift;

  my $exists = 0;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT xref_id FROM primary_xref WHERE xref_id=?");
  $sth->execute($xref_id) || die $dbi->errstr;
  my @row = $sth->fetchrow_array();
  my $result = $row[0];
  $exists = 1 if (defined $result);

  return $exists;

}

# --------------------------------------------------------------------------------

1;

package BaseParser;

use strict;

use DBI;
use Digest::MD5 qw(md5_hex);
use POSIX qw(strftime);

use SwissProtParser;

my $host = "ecs1g";
my $port = 3306;
my $database = "glenn_test_xref";
my $user = "ensadmin";
my $password = "ensembl";

my $base_dir = ".";

my $dbi;
my %dependent_sources;
my %taxonomy2species_id;

my %filetype2parser = (
		       "UniProtSwissProt" => "SwissProtParser",
		       "UniProtTrEMBL"    => "SwissProtParser"
		      );

run() if (!defined(caller()));

# --------------------------------------------------------------------------------
# Get info about files to be parsed from the database

sub run {

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT * FROM source WHERE url IS NOT NULL");
  $sth->execute();
  my ($source_id, $name, $url, $checksum, $modified_date, $upload_date, $release);
  $sth->bind_columns(\$source_id, \$name, \$url, \$checksum, \$modified_date, \$upload_date, \$release);
  while (my @row = $sth->fetchrow_array()) {

    # Download each source into the appropriate directory for parsing later
    # Convention is that part of name up to _ is file type
    my ($type) = $name =~ /^([^_]+)_.*/;
    my $dir = $base_dir . "/" . $type;
    mkdir $dir if (!-e $dir);

    my ($file) = $url =~ /.*\/(.*)/;

    print "Downloading $url to $dir/$file\n";
    # TODO remove --timestamping??
    my $result = system("wget", "--quiet", "--timestamping", "--directory-prefix=$dir", $url);

    # compare checksums and parse/upload if necessary
    # need to check file size as some .SPC files can be of zero length
    my $file_cs = md5sum("$dir/$file");
    if (!defined $checksum || $checksum ne $file_cs) {

      if (-s "$dir/$file") {

	print "Checksum for $file does not match, parsing\n";

	#update_source($dbi, $source_id, $file_cs, $file);

	my $parserType = $filetype2parser{$type};
	print "Parsing $file with $parserType\n";
	#$parserType->run("$dir/$file", $source_id);

      } else {
	
	print $file . " has zero length, skipping\n";

      }

    } else {

      print $file . " has not changed, skipping\n";

    }
  }

}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "BaseParser";

  return $self;

}

# --------------------------------------------------------------------------------
# Get source ID for a particular source name

sub get_source_id {

  my ($self, $name) = @_;

  my $sth = $dbi->prepare("SELECT source_id FROM source WHERE name='" . $name . "'");
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $source_id;
  if (defined @row) {
    $source_id = $row[0];
  } else {
    warn("Couldn't get source ID for name $name\n");
    $source_id = -1;
  }

  return $source_id;

}

# --------------------------------------------------------------------------------
# Upload xrefs to the database

sub upload_xrefs {

  my ($self, @xrefs) = @_;

  my $dbi = dbi();

  # remove all existing xrefs with same source ID
  my $source_id = $xrefs[0]->{SOURCE_ID};
  my $del_sth = $dbi->prepare("DELETE FROM xref WHERE source_id=$source_id");
  $del_sth->execute();

  # upload new ones
  my $xref_sth = $dbi->prepare("INSERT INTO xref (accession,label,source_id,species_id) VALUES(?,?,?,?)");
  my $pri_sth = $dbi->prepare("INSERT INTO primary_xref VALUES(?,?,?,?,?)");
  my $syn_sth = $dbi->prepare("INSERT INTO synonym VALUES(?,?,?)");
  my $dep_sth = $dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");

  foreach my $xref (@xrefs) {

    # Create entry in xref table and note ID
    $xref_sth->execute($xref->{ACCESSION},
		       $xref->{LABEL},
		       $xref->{SOURCE_ID},
		       $xref->{SPECIES_ID}) || die $dbi->errstr;

    # get ID of xref just inserted
    my $xref_id = $xref_sth->{'mysql_insertid'};

    # create entry in primary_xref table with sequence
    # TODO experimental/predicted????
    $pri_sth->execute($xref_id,
		      $xref->{SEQUENCE},
		      'peptide',
		      'experimental',
		      $xref->{SOURCE_ID}) || die $dbi->errstr;

    # if there are synonyms, create xrefs for them and entries in the synonym table
    foreach my $syn (@{$xref->{SYNONYMS}}) {

      $xref_sth->execute($syn,
			 $xref->{LABEL},
			 $xref->{SOURCE_ID},
			 $xref->{SPECIES_ID}) || die $dbi->errstr;

      my $syn_xref_id = $xref_sth->{'mysql_insertid'};
      $syn_sth->execute($xref_id, $syn_xref_id, $xref->{SOURCE_ID} ) || die $dbi->errstr;

    } # foreach syn

    # if there are dependent xrefs, add xrefs and dependent xrefs for them
    foreach my $depref (@{$xref->{DEPENDENT_XREFS}}) {

      my %dep = %$depref;

      $xref_sth->execute($dep{ACCESSION},
			 $xref->{LABEL},
			 $xref->{SOURCE_ID},
			 $xref->{SPECIES_ID}) || die $dbi->errstr;

      my $dep_xref_id = $xref_sth->{'mysql_insertid'};
			
      $dep_sth->execute($xref_id, $dep_xref_id, '', $dep{SOURCE_ID} ) || die $dbi->errstr;
      # TODO linkage anntation?

    } # foreach dep


  } # foreach xref

  $del_sth->finish() if defined $del_sth;
  $xref_sth->finish() if defined $xref_sth;
  $pri_sth->finish() if defined $pri_sth;

}

# --------------------------------------------------------------------------------
# Get & cache a hash of all the source names for dependent xrefs (those that are
# in the source table but don't have an associated URL etc)

sub get_dependent_xref_sources {

  my $self = shift;

  if (!defined %dependent_sources) {

    my $dbi = dbi();
    my $sth = $dbi->prepare("SELECT name,source_id FROM source WHERE url IS NULL");
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
# Update a row in the source table

sub update_source {

  my ($dbi, $source_id, $checksum, $file) = @_;
  open(FILE, $file);
  my $file_date = POSIX::strftime('%Y%m%d%H%M%S', localtime((stat($file))[9]));
  close(FILE);

  my $sql = "UPDATE source SET checksum='" . $checksum . "', file_modified_date='" . $file_date . "', upload_date=NOW() WHERE source_id=" . $source_id;
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

1;

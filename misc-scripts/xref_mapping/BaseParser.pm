package BaseParser;

use DBI;

my $host = "ecs1g";
my $port = 3306;
my $database = "glenn_test_xref";
my $user = "ensadmin";
my $password = "ensembl";

my $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$database",
		       "$user",
		       "$password",
		       {'RaiseError' => 1}) || die "Can't connect to database";

# --------------------------------------------------------------------------------
# Upload a source object to the source table.
# Return the ID of the new row in the database.

sub upload_source {

  my ($self, $source) = @_;

  my $sth = $dbi->prepare("INSERT INTO source (name,url,file_modified_date,upload_date,release) VALUES(?,?,?, NOW(),?)");
  $sth->execute($source->{NAME},
		$source->{URL},
		$source->{FILE_MODIFIED_DATE}, "") || die $dbi->errstr;

  $sth = $dbi->prepare("SELECT last_insert_id()");
  $sth->execute();
  my $id;
  while(my @row = $sth->fetchrow_array()) {
    $id = $row[0];
  }

  $sth->finish() if defined $sth;

  return $id;

}

# --------------------------------------------------------------------------------
# Upload xrefs to the database

sub upload_xrefs {

  my ($self, @xrefs) = @_;

  my $sth;
  foreach my $xref (@xrefs) {

    # Create entry in xref table and note ID
    $sth = $dbi->prepare("INSERT INTO xref (accession,label,source_id,species_id) VALUES(?,?,?,?)");
    $sth->execute($xref->{ACCESSION},
		  $xref->{LABEL},
		  $xref->{SOURCE_ID},
		  $xref->{SPECIES_ID}) || die $dbi->errstr;

    # get ID of xref just inserted
    $sth = $dbi->prepare("SELECT last_insert_id()");
    $sth->execute();
    my $id;
    while(my @row = $sth->fetchrow_array()) {
      $id = $row[0];
    }

    # create entry in primary_xref table with sequence
    # TODO experimental/predicted????
    $sth = $dbi->prepare("INSERT INTO primary_xref VALUES(?,?,?,?,?)");
    $sth->execute($id,
		  $xref->{SEQUENCE},
		  'peptide',
		  'experimental',
		  $xref->{SOURCE_ID}) || die $dbi->errstr;

  }

  $sth->finish() if defined $sth;

}

# --------------------------------------------------------------------------------

sub dbi {

  my $self = shift;

  return $dbi;

}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "BaseParser";
  return $self;

}

# --------------------------------------------------------------------------------

1;

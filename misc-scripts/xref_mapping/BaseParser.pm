package BaseParser;

use DBI;
use strict;

my $host = "ecs1g";
my $port = 3306;
my $database = "glenn_test_xref";
my $user = "ensadmin";
my $password = "ensembl";

my $dbi;
my %dependent_sources;

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

  my $self= shift;

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


# --------------------------------------------------------------------------------

1;

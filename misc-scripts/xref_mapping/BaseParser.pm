package BaseParser;

use DBI;
use strict;

my $host = "ecs1g";
my $port = 3306;
my $database = "glenn_test_xref";
my $user = "ensadmin";
my $password = "ensembl";

my $dbi;

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

  # remove all existing xrefs with same source ID
  my $source_id = $xrefs[0]->{SOURCE_ID};
  my $del_sth = $dbi->prepare("DELETE FROM xref WHERE source_id=$source_id");
  $del_sth->execute();

  # upload new ones
  my $xref_sth = $dbi->prepare("INSERT INTO xref (accession,label,source_id,species_id) VALUES(?,?,?,?)");
  my $pri_sth = $dbi->prepare("INSERT INTO primary_xref VALUES(?,?,?,?,?)");
  my $syn_sth = $dbi->prepare("INSERT INTO synonym VALUES(?,?,?)");

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

  } # foreach xref

  $del_sth->finish() if defined $del_sth;
  $xref_sth->finish() if defined $xref_sth;
  $pri_sth->finish() if defined $pri_sth;

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

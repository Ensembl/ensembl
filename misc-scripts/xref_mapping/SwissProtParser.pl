# Parse SwissProt files to create xrefs.

use strict;
use DBI;
use Data::Dumper;

my $host = "ecs1g";
my $port = 3306;
my $database = "glenn_test_xref";
my $user = "ensro";
my $password = "";

# --------------------------------------------------------------------------------
# Parse command line

if (scalar(@ARGV) != 1) {
  print "\nUsage: SwissProtParser.pl file.SPC\n\n";
  exit(1);
}

my $file = $ARGV[0];

# --------------------------------------------------------------------------------
# Parse file into array of xref objects

open(SWISSPROT, $file) || die "Can't open Swissprot file $file\n"; 

my @xrefs;

my $previous_du = $/;
$/ = "\/\/\n";

while (<SWISSPROT>) {

  my $xref;
  ($xref->{ACCESSION}) =$_ =~ /AC\s+(\w+);/;
  ($xref->{LABEL}) = $_ =~ /DE\s+(.+)/;

  #print $xref->{ACCESSION} . " " . $xref->{LABEL} ."\n";

  push @xrefs, $xref;

}

$/ = $previous_du;

print "Read " . scalar(@xrefs) ." xrefs from $file\n";

# --------------------------------------------------------------------------------
# Create source object to be loaded into source table

my $source;
$source = {
	   NAME => "SwissProt",
	   URL  => $file,
	   UPLOAD_DATE => time(),
	   FILE_MODIFIED_DATE => (stat($file))[9]
};

#print Dumper($source);

# TODO - dates as formatted strings
# TODO URL? Release?

# --------------------------------------------------------------------------------

my $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$database", "$user", "$password", {'RaiseError' => 1}) || die "Can't connect to database";

# upload the source and get the ID
#my $sth = $dbi->prepare("INSERT INTO source (name,url,file_modified_date,upload_date,release) VALUES(?,?,?,?,?,?)");
#$sth->execute($source->{NAME}, $source->{URL}, $source->{FILE_MODIFIED_DATE}, $source->{UPLOAD_DATE}, "") || die $dbi->errstr;
# TODO last_insert_id() in mySQL
# TODO better error handling

# TODO what about old versions?

# TODO upload xrefs, source
# TODO set source_id in each xref appropriately

# Dump primary xrefs to FASTA file

use strict;
use DBI;

my $file = "xref_dna.fasta";

my $host = "ecs1g";
my $port = 3306;
my $database = "glenn_test_xref";
my $user = "ensadmin";
my $password = "ensembl";

my $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$database",
		       "$user",
		       "$password",
		       {'RaiseError' => 1}) || die "Can't connect to database";

open(FILE, ">" . $file);

my $sth = $dbi->prepare("SELECT x.xref_id, px.sequence FROM primary_xref px, source so, xref x, species sp WHERE sp.name='homo_sapiens' AND so.name='RefSeq' AND x.species_id=sp.species_id AND px.source_id=so.source_id AND x.xref_id=px.xref_id AND px.sequence_type='dna' LIMIT 500");
$sth->execute();

my ($xref_id, $sequence);
$sth->bind_columns(\$xref_id, \$sequence);

while (my @row = $sth->fetchrow_array()) {

  print FILE ">$xref_id\n$sequence\n";

}

close(FILE);


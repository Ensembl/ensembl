use strict;
use warnings;
use IO::File;
use DBI;

my $input = 'QTL2.txt';
my $host = '127.0.0.1';
my $dbname = 'rattus_norvegicus_core_13_2';
my $port   = '3308';
my $user   = 'ensro';
my $pass   = '';

my $db = DBI->connect("DBI:mysql:host=$host;dbname=$dbname;port=$port;", $user, $pass, {'RaiseError' => 1});

my $q = qq(
SELECT source_primary_id, qtl_id FROM qtl WHERE source_database = 'rat genome database'
);
my $sth = $db->prepare($q);
$sth->execute();

my %rgd_hash = map {@$_} @{$sth->fetchall_arrayref};

$sth->finish();
$db->disconnect();

my $fh = new IO::File;
$fh->open($input);

<$fh>; #throw away heading row
my %ratmap_hash;
while(<$fh>) {
  chomp;
  my @a = split(/\t/, $_);
  my ($rgd_id, $ratmap_id) = ($a[4], $a[3]);
  
  next if(!$ratmap_id || !$rgd_id);

  if(exists($ratmap_hash{$ratmap_id})) {
    warn "DUPLICATE identifier $ratmap_id\n";
    next;
  }

  $ratmap_hash{$ratmap_id} = 1;

  if(!exists($rgd_hash{$rgd_id})) {
    warn "QTL RGD:$rgd_id not in database\n";
    next;
  }

  my $ens_id = $rgd_hash{$rgd_id};
  if($ens_id) {
    print "INSERT INTO qtl_synonym (qtl_id, source_database, source_primary_id) VALUES ($ens_id, 'ratmap', $ratmap_id);\n";
  } else {
    #print "#RatMap QTL $ratmap_id not found\n";
  }
}

$fh->close();

foreach my $rgd_id (keys %rgd_hash) {
  my $ens_id = $rgd_hash{$rgd_id};
  print "INSERT INTO qtl_synonym (qtl_id, source_database, source_primary_id) VALUES ($ens_id, 'rat genome database', $rgd_id);\n";
}


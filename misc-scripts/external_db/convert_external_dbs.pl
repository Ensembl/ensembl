use strict;
use warnings;
use DBI;

my $user = 'ecs2dadmin';
my $host = 'ecs2d';
my $pass = 'TyhRv';
my $dbname = shift;

my $dbh = DBI->connect("DBI:mysql:host=$host;dbname=$dbname;", $user, $pass,
		       {RaiseError => 1});

#
# Store the old external_db table in a hash
#
my $sth = $dbh->prepare('SELECT external_db_id, db_name
                         FROM external_db');

print STDERR "READING OLD EXTERNAL DB\n";
$sth->execute();
my %old_ext_db = map {$_->[0], $_->[1]} @{$sth->fetchall_arrayref};
$sth->finish();


#
# drop the existing external_db table and replace it with the new table
#
print STDERR "REPLACING OLD EXTERNAL DB\n";

`cat external_db.sql | mysql -h $host -u $user -p$pass $dbname`;

#
# Store the new external_db table in a hash
#
print STDERR "READING NEW EXTERNAL DB\n";
$sth->execute();
my %new_ext_db = map {$_->[1], $_->[0]} @{$sth->fetchall_arrayref};
$sth->finish();

#
# update each row in the xref table
#
print STDERR "UPDATING XREF TABLE\n";
$sth = $dbh->prepare('SELECT external_db_id, xref_id FROM xref');
my($external_db_id, $xref_id);
$sth->execute();
$sth->bind_columns(\$external_db_id, \$xref_id);

my $update_sth = 
  $dbh->prepare('UPDATE xref SET external_db_id = ? WHERE xref_id = ?');

my $count = 0;

while($sth->fetch()) {
  my $dbname = $old_ext_db{$external_db_id};
  my $id     = $new_ext_db{$dbname};
  if($id) {
    $update_sth->execute($id, $xref_id);
    $count++;
    if($count % 1_000 == 0) {
      print STDERR '.';
    }
  } else {
    warn("Could not convert ext_id=[$external_db_id] dbname=[$dbname]\n");
  }
}

$sth->finish();

$dbh->disconnect();

print STDERR "COMPLETE.  Converted $count xrefs.\n";


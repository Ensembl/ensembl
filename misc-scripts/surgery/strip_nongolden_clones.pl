#!/usr/local/bin/perl


use DBI;
use strict;

my $host   = 'ecs1e';
my $dbname = 'ens_apr01_gb';
my $dbuser = 'ensro';
my $password = undef;

my $dsn = "DBI:mysql:database=$dbname;host=$host;";

my $dbh = DBI->connect("$dsn","$dbuser",$password, {RaiseError => 1});


# select internal ids for all the clones

my $sth = $dbh->prepare("select internal_id from clone;");
$sth->execute;

my %clone;
my $count = 0;
while( my ($internal_id) = $sth->fetchrow_array ) {
  $clone{$internal_id} = 0;
  $count++;
}
print STDERR "Got $count clones\n";


# select internal ids for all clones joined to the static golden path

my $sth = $dbh->prepare("select distinct(c.clone) from contig c,static_golden_path s where c.internal_id = s.raw_id;");
$sth->execute;

my $gcount = 0;
while( my ($internal_id) = $sth->fetchrow_array ) {
  $clone{$internal_id} = 1;
  $gcount++;
}

print STDERR "Got $gcount golden clones\n";



# foreach clone, if 0, find contig and dna identifiers and issue delete

foreach my $clone_id ( keys %clone ) {
  if( $clone{$clone_id} == 1 ) {
    next;
  }
  
  print "delete from clone where internal_id = $clone_id;\n";
  
  my $sth = $dbh->prepare("select dna,internal_id from contig where clone = $clone_id");
  $sth->execute;

  while( my ($dna,$contig) = $sth->fetchrow_array ) {
    print "delete from dna where id = $dna;\n";
    print "delete from contig where internal_id = $contig;\n";
  }

}



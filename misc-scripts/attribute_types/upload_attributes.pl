#!/usr/local/ensembl/bin/perl -w
#
# Stolen from external db directory this script performs
# the similar function of uploading all valid attribute types
# into the attrib_type table of all databases that are going to
# be released.

use strict;

use Getopt::Long;
use DBI;
use IO::File;
use FindBin;

my (
  $host,        $user,     $pass,
  $port,        @dbnames,  $file,
  $release_num, $nobackup, $forceTableWrite
);

GetOptions(
  "dbhost|host=s",     \$host,
  "dbuser|user=s",     \$user,
  "dbpass|pass=s",     \$pass,
  "dbport|port=i",     \$port,
  "file=s",            \$file,
  "nobackup",          \$nobackup,
  "dbname|dbnames=s@", \@dbnames,
  "release_num=i",     \$release_num,
  "force_table_write", \$forceTableWrite
);

$file ||= $FindBin::Bin . "/attrib_type.txt";
usage() if ( !$host );
usage()
  if ( ( $release_num && @dbnames ) || ( !$release_num && !@dbnames ) );

#release num XOR dbname are required
$port ||= 3306;

my $dsn = "DBI:mysql:host=$host;port=$port";
my $db = DBI->connect( $dsn, $user, $pass, { RaiseError => 1 } );

if ($release_num) {
  @dbnames =
    map { $_->[0] } @{ $db->selectall_arrayref("show databases") };

  #
  # filter out all non-core databases
  #
  @dbnames = grep {
/^[a-zA-Z]+\_[a-zA-Z]+\_(core|otherfeatures|cdna|vega)\_${release_num}\_\d+[A-Za-z]?$/
  } @dbnames;
}

#
# make sure the user wishes to continue
#
print STDERR "The following databases "
  . "will have their attrib_type tables updated:\n  ";
print join( "\n  ", @dbnames );
print "\ncontinue with update (yes/no)>  ";

my $input = lc(<STDIN>);
chomp($input);
if ( $input ne 'yes' ) {
  print "Attrib_type loading aborted\n";
  exit();
}

my $attribs = read_attrib_file($file);

my $table_consistent;
# if any attrib_types are loaded that are different from the file, a
# consistency problem is reported and the upload is not done.
for my $database (@dbnames) {

  if ( !$nobackup ) {
    backup_attribute_types( $host, $user, $pass, $port, $database );
  }

  $table_consistent = check_consistency( $attribs, $database, $db );

  # This has been introduced in e54: when we are sure that the data in
  # the file is consistent (e.g. some attributes have been removed from
  # the file since they are not longer are needed) you should force the
  # writing option.

  if ( $table_consistent || $forceTableWrite ) {
    # consistent
    $db->do("use $database");
    $db->do("delete from attrib_type");

    load_attribs( $db, $attribs );
  } else {
    print STDERR "Repairing $database, not consistent!\n";
    repair( $attribs, $database, $db );
    print STDERR "If you are sure the file is up to date, "
      . "you can use the --force_table_write"
      . "to overwrite the information in the attrib_table, "
      . "check RelCoord for more information\n";
  }
} ## end for my $database (@dbnames)

# Move attrib types wih the same code to the common attrib_type table
# ones that are not in the table move to an attrib_type_id that is not
# used.
sub repair {

  my ( $attribs, $database, $db ) = @_;

  $db->do("use $database");

  my @tables = qw( seq_region_attrib misc_attrib translation_attrib
    transcript_attrib gene_attrib);
  my $ref = $db->selectall_arrayref("show create table attrib_type");
  my $create_table = $ref->[0]->[1];

  $db->do("alter table attrib_type rename old_attrib_type");
  $db->do($create_table);

  load_attribs( $db, $attribs );

  $db->do( "delete oat "
      . "from old_attrib_type oat, attrib_type at "
      . "where oat.attrib_type_id = at.attrib_type_id "
      . "and oat.code = at.code" );

  # what remains in old attrib type ?
  #  Entries with a code that is unknown in general file and that
  #  shouldn't really happen. If it happens, the code needs to be
  #  appended to attrib_type table and the attrib type_ids will be
  #  updated in the feature tables.

  #  Entries with a code that is known, but has different
  #  attrib_type_id. Feature tables will be updated.

  $db->do( "create table tmp_attrib_types "
      . "select oat.attrib_type_id, oat.code, oat.name, oat.description "
      . "from old_attrib_type oat "
      . "left join attrib_type at "
      . "on oat.code = at.code "
      . "where at.code is null" );
  $db->do( "insert into attrib_type( code, name, description) "
      . "select code, name, description "
      . "from tmp_attrib_types" );

  $ref = $db->selectall_arrayref("select code from tmp_attrib_types");
  $db->do("drop table tmp_attrib_types");

  print STDERR "Database $database todo:\n";
  if (@$ref) {
    print STDERR "  Missing codes ",
      join( ", ", map { $_->[0] } @$ref ), "\n";
  }

  my %missing_codes = map { $_->[0], 1 } @$ref;

  $ref =
    $db->selectall_arrayref("select code from old_attrib_type oat ");

  my @updated_codes;
  for my $code_ref (@$ref) {
    if ( !exists $missing_codes{ $code_ref->[0] } ) {
      push( @updated_codes, $code_ref->[0] );
    }
  }

  print STDERR "  Updated codes ", join( ", ", @updated_codes ), "\n";

  # now do multi table updates on all tables
  for my $up_table (@tables) {
    $db->do( "update $up_table tb, attrib_type at, old_attrib_type oat "
        . "set tb.attrib_type_id = at.attrib_type_id "
        . "where tb.attrib_type_id = oat.attrib_type_id "
        . "and oat.code = at.code " );
  }

  $db->do("drop table old_attrib_type");
} ## end sub repair

sub load_attribs {
  my ( $db, $attribs ) = @_;
  my $sth;
  $sth = $db->prepare(
    "insert into attrib_type( attrib_type_id, code, name, description) "
      . "values(?,?,?,?)" );
  for my $attrib (@$attribs) {
    $sth->execute(
      $attrib->{'attrib_type_id'}, $attrib->{'code'},
      $attrib->{'name'},           $attrib->{'description'} );
  }
}

# alternatively consistency can be enforced to a certain degree
sub check_consistency {
  my $attribs  = shift;
  my $database = shift;
  my $db       = shift;

  my ( %db_codes, %file_codes );
  map { $file_codes{ $_->{'attrib_type_id'} } = $_->{'code'} }
    @$attribs;

  $db->do("use $database");
  my $sth =
    $db->prepare( "SELECT attrib_type_id, code, name, description "
      . "FROM attrib_type" );
  $sth->execute();
  while ( my $arr = $sth->fetchrow_arrayref() ) {
    $db_codes{ $arr->[0] } = $arr->[1];
  }

  # check if any ids in the database collide with the file
  my $consistent = 1;
  for my $dbid ( keys %db_codes ) {
    if ( !exists $file_codes{$dbid} ) {
      printf( "Not consistent: code '%d' ('%s') not in file.\n",
        $dbid, $db_codes{$dbid} );
      $consistent = 0;
    } elsif ( $file_codes{$dbid} ne $db_codes{$dbid} ) {
      printf( "Not consistent: code '%d' is '%s' in file "
          . "but '%s' in database.\n",
        $dbid, $file_codes{$dbid}, $db_codes{$dbid} );
      $consistent = 0;
    }
  }

  return $consistent;
} ## end sub check_consistency

sub read_attrib_file {
  my $file = shift;
  #
  # read all attrib_type entries from the file
  #
  my $fh = IO::File->new();
  $fh->open($file) or die("could not open input file $file");

  my @rows;
  my $row;
  while ( $row = <$fh> ) {
    chomp($row);
    next if ( $row =~ /^\S*$/ );
    next if ( $row =~ /^\#/ );

    my @a = split( /\t/, $row );

    push @rows,
      {
      'attrib_type_id' => $a[0],
      'code'           => $a[1],
      'name'           => $a[2],
      'description'    => $a[3] };
  }
  $fh->close();
  return \@rows;
} ## end sub read_attrib_file

sub backup_attribute_types {
  my ( $host, $user, $pass, $port, $dbname ) = @_;

  unless (
    system(
          "mysql -h$host -P$port -u$user -p'$pass' -N "
        . "-e 'select * from attrib_type' $dbname "
        . "> $dbname.attrib_type.backup; "
        . "gzip -9 -f $dbname.attrib_type.backup"
    ) == 0
    )
  {
    print STDERR "Can't dump the original attrib_type table "
      . "from $dbname for backup\n";
    exit 1;
  } else {
    print STDERR "Original attrib_type table backed up in "
      . "$dbname.attrib_type.backup\n";
  }
}

sub usage {
  GetOptions(
    "host=s",            \$host,
    "user=s",            \$user,
    "pass=s",            \$pass,
    "port=i",            \$port,
    "file=s",            \$file,
    "dbnames=s@",        \@dbnames,
    "release_num=i",     \$release_num,
    "nobackup",          \$nobackup,
    "force_table_write", \$forceTableWrite
  );

  print STDERR <<EOC;
Usage:
  perl upload_attributes.pl --host ... --user ... --port ... --pass ...

  and either

    --dbnames homo_sap_1 --dbname homo_sap_2 --dbnames ...

  or

   --release_num ...

  Use --file ... if not using the default 'attrib_type.txt'.

  The new option force_table_write has been introduced: when you run
  the script and there is inconsistency between the database and the
  file (there is the "not consistent!!" in the output), you should first
  ensure if the consistency is expected (some attributes were removed
  from the file because they are not longer needed), or not expected
  (someone introduced a new attribute_type in the database but forgot to
  update the attrib_type.txt file). In the first situation, you should
  use the --force_table_write when running the script. In the second
  situation, the new attrib_type must be added to the file and run it
  again without the force flag

EOC
  exit;
} ## end sub usage

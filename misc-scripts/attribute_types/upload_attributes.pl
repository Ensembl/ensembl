# Stolen from external db directory this script performs
# the similar function of uploading all valid attribute types
# into the attrib_type table of all databases that are going to
# be released.


use strict;

use Getopt::Long;
use DBI;
use IO::File;
use FindBin;


my ( $host, $user, $pass, $port,@dbnames, $file, $release_num);

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "file=s", \$file,
            "dbnames=s@", \@dbnames,
	    "release_num=i", \$release_num
	  );

$file ||= $FindBin::Bin."/attrib_type.txt";
usage() if(!$host);



#release num XOR dbname are required
usage() if(($release_num && @dbnames) || (!$release_num && !@dbnames));
my $attribs = read_attrib_file( $file );

$port ||= 3306;

my $dsn = "DBI:mysql:host=$host;port=$port";

my $db = DBI->connect( $dsn, $user, $pass, {RaiseError => 1} );

if($release_num) {
  @dbnames = map {$_->[0] } @{ $db->selectall_arrayref( "show databases" ) };
  
  #
  # filter out all non-core databases
  #
  @dbnames = grep {/^[a-zA-Z]+\_[a-zA-Z]+\_(core|est|estgene|vega)\_${release_num}\_\d+[A-Za-z]?$/} @dbnames;
}


#
# make sure the user wishes to continue
#
print STDERR "The following databases will be external_db updated:\n  ";
print join("\n  ", @dbnames);
print "\ncontinue with update (yes/no)>  ";

my $input = lc(<STDIN>);
chomp($input);
if($input ne 'yes') {
  print "Attrib_type loading aborted\n";
  exit();
}

# if any attrib_types are loaded that are different from 
# the file, a consistency problem is reported and the 
# upload is not done.

# alternatively consistency can be enforceed to a certain degree 
sub check_consistency {
  my $attribs = shift;
  my $database_names = shift;
  my $db = shift;

}



sub read_attrib_file {
  my $file = shift;
  #
  # read all attrib_type entries from the file
  #
  my $fh = IO::File->new();
  $fh->open($file) or die("could not open input file $file");
  
  my @rows;
  my $row;
  while($row = <$fh>) {
    chomp($row);
    next if( /^\S*$/ );
    next if ( /^\#/ );

    my @a = split(/\t/, $row);
    
    push @rows, {'attrib_type_id' => $a[0],
		 'code' => $a[1],
		 'name' => $a[2],
		 'description'  => $a[3]};
  }
  $fh->close();
  return \@rows;
}



# loading is:
# which clones to load read into hash
# generate a table for putting into mysql

# check the next free number in dna table.

# generates flat files for upload into db
# needs thorough editing on filennames before usage
# input is a clonelist file and *ffa* files from the freeze

# look at the upload subroutine as well. just three tables need
# uploading


use DBI;

  

$dsn = "dbi:mysql:database=ensembl_freeze17;host=localhost";
$db = DBI->connect( $dsn, 'root' );

if( shift eq '-upload' ) {
  upload();
  exit;
}

$sth = $db->prepare( "select max(id) from dna" );
$sth->execute;
$maxId = ($sth->fetchrow_array)[0];

$sth = $db->prepare( "select now()" );
$sth->execute;
$now = ($sth->fetchrow_array)[0];

my $nextDNAid = $maxId+1;

my %clonehash;

$sourceDir = '/nfs/disk100/humpub/th/unfinished_ana/ncbi_files/freeze_jul17';
@files = glob( "$sourceDir/*ffa*" );
open( DNA, ">dna.txt" );
open( CONTIG, ">contig.txt" );
open( CLONE, ">clone.txt" );
open( LOAD, "clonesToLoad.txt" );

while( <LOAD> ) {
  chomp;
  $clonehash{$_} = 1;
}

print ( "Clones to load: ",scalar( keys %clonehash ), "\n" );
$clonecount = 0;

for $file (@files) {
  open( FH, "gzip -d -c $file |" ) or 
    dir( "couldnt read $file" );
  $load = 0;
  while( <FH> ) {
    if( /^>(\S*) \(([^\.]*)\.(\d+):(\d+)\.\.(\d+)\)/ ) {
      if( $load ) {
        print DNA "$dna\t$now\t$nextDNAid\n";
	$nextDNAid++;
      }
      $dna = '';

      if( $clonehash{$2} == 1 ) {
        # clone is not yet written to file
	print CLONE "$2\t\t$3\t$now\t$now\t$now\n";
	$clonecount++;
	if( $clonecount % 100 == 0 ) {
	  print( "$clonecount clones loaded.\n" );
	}
      }

      if( defined $clonehash{$2} ) {
        $load = 1;
	my ( $international_id, $cloneId, $offset, $end ) = ( $1, $2, $4, $5 ); 
	my $length = $end - $offset + 1; 
	
	my $num = $clonehash{$cloneId};
	$clonehash{$cloneId}++;
	my $id = sprintf( "%s.%05d", $cloneId, $num );
	print CONTIG "$id\t$cloneId\t$length\t$offset\t$nextDNAid\t0\t$international_id\n";
      } else {
        $load = 0;
      }
      
    } else {
      if( $load ) {
        if( /([\S]*)/ ) {
          $dna .= uc( $1 );
        }  
      }
    }
  }
  close FH;
  if( $load ) {
    print DNA "$dna\t$now\t$nextDNAid\n";
    $nextDNAid++;
  }
}

sub upload {
# load stuff
  my $rv;
  
#  print "loading clones\n";
#  $rv = $db->do( "load data infile '/mysql/michele/clone.txt' into table clone( id, embl_id, embl_version, created, modified, stored )" );
#  print "Result: $rv ", $DBI::errstr ,"\n";
#  print "loading contigs\n";
#  $rv = $db->do( "load data infile '/mysql/michele/contig.txt' into table contig( id, clone, length, offset, dna, chromosomeId, international_id )" );
#  print "Result: $rv ", $DBI::errstr ,"\n";
  print "loading dna\n";
  open( DNA, '/mysql/michele/dna.txt' ) or die "couldnt open dna file";
  my $count = 0;
 
  while( $line = <DNA> ) {
    my @arr = split( /[\t]/, $line );
    $rv = $db->do( "insert into dna set sequence='$arr[0]',created='$arr[1]',id=$arr[2]" );    
#   print( "insert into dna set sequence='$arr[0]',created='$arr[1]',id=$arr[2] \n" );    
#  $rv = $db->do( "load data infile '/mysql/michele/dna.txt' into table dna( sequence, created, id )" );
    if( $rv != 1 ) {
      print "Result: $rv ", $DBI::errstr ,"\n";
    }
    $count++;
    if( $count % 500 == 0 ) {
      print "$count dna loaded\n";
    }
  } 
}

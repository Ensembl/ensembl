# sql to do
use DBI;

# give as argument a file with clonelist.
# It deletes that clones from the database
# provide database details further down in the script 
# it deletes all predictions and overlaps
# maybe you have to add more delete from lines ...
# definatly feature table needs some deletion, I just wasnt bothering

# read in the clone ids to delete
open( DEL, shift ) or die ( "couldnt open clonelist" );
my @clonelist;
while( <DEL> ) {
  chomp;
  push( @clonelist, $_ );
}

$inClause1 = "('".join( "','", @clonelist )."')";

$dsn = "dbi:mysql:database=ensembl_freeze17;host=localhost";
$db = DBI->connect( $dsn, 'root' );

eval {
$db->do( "delete from exon" );
$db->do( "delete from exon_pair");
$db->do( "delete from exon_transcript" );
$db->do( "delete from contigoverlap" );
};

$rv = $db->do( "delete from clone where id in $inClause1" );
print( "delete clone returned $rv\n" );
print( "Clonelist contained ", scalar( @clonelist), " elements.\n" );
print( "Error (if any): ",$DBI::errstr, "\n" );

# print ( "delete from clone where id in $inClause1" );

$sth1 = $db->prepare( "select dna from contig where clone in $inClause1" );
$rv = $sth1->execute;
print( "Error (if any): ",$DBI::errstr, "\n" );
my @dnalist;

while( @arr = $sth1->fetchrow_array ) {
  push( @dnalist, $arr[0] );
}

$inClause2 = "(".join( ",", @dnalist ).")";
 
# $rv = $db->do( "delete from dna where id in $inClause2" );
print( "delete from dna $rv\n" );
print( "Error (if any): ",$DBI::errstr, "\n" );
$rv = $db->do( "delete from contig where clone in $inClause1" );
print( "delete from contig $rv\n" );
print( "Error (if any): ",$DBI::errstr, "\n" );

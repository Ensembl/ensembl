
# This script takes as argument the name of a gzipped cigar file, 
# which contains exonerate result from mouse traces to human dna

# you have to provide an empty ensembl database with just contig table 
# filled 
# lots of empty diskspace
# change db deatils below !!



use DBI;

# only use with local host!
my $file = shift;

my $dsn = "DBI:mysql:database=mouse_jun_hum_apr_arne;host=localhost";
my $dbh = DBI->connect($dsn,'root');
my $tmpfile = "/scratch1/ensembl/arne/mouse_features.txt";




open( EXONERATE, "gzip -dc $file|" ) or die( "cant open exonerate cigar line file" );
open( TMP, ">$tmpfile" ) or die ("cant open tmpfile $tmpfile");

# print join( "\t", contig, seq_start, seq_end, score, strand, analysis, "exonerate_gapped", hstart, hend, hid )

$sth = $dbh->prepare( "select internal_id, id from contig" );

my ( $internal, $external );
$sth->execute();
$sth->bind_columns( \$internal, \$external );

while( $sth->fetch() ) {
  $contigHash{$external} = $internal;
}

my $count = 0;

while( <EXONERATE> ) {
  chomp;
  split;
  ($hid) = $_[2] =~ /(.*)-[0-9]*$/;
  $hstart = $_[3]+1;
  $hend = $_[4];
  $strand = ($_[5] eq '-') ? -1 : 1;
  $contig = $contigHash{$_[6]};
  $seq_start = $_[7]+1;
  $seq_end = $_[8];
  $score = $_[10];
  $analysis = 1;

  print TMP ( join( "\t",( $contig, $seq_start, $seq_end, $score, 
			   $strand, $analysis, "exonerate_gapped", 
			   $hstart, $hend, $hid )), "\n" );
  $count++;
  if( $count % 100000 == 0 ) {
    print STDERR ".";
#    exit;
  }
}


$dbh->do( "load data infile '$tmpfile' into table feature( contig, seq_start, seq_end, score, strand, analysis, name, hstart, hend, hid )" );
$dbh->do( "insert into analysis( id, db, db_version, program, program_version, gff_source, gff_feature) values( 1, \"exonerate_gapped\", \"\", \"exonerate\", \"1\", \"exonerate_gapped\", \"similarity\" )" );
$dbh->disconnect();

exit;

  
# algorithm		     

# build contig internal hash
# read cigar lines from gzip -dc pipe
# print feature lines into tmpfile

# issue load data infile command .
# fill in one row in analysis table..


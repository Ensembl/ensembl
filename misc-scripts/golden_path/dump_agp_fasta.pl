# load the ncbi ffa header lines
# read the clone versions
# update the versions in the db
# one argument. the agp filename
# db connection is hardcoded


use DBI;
use Bio::SeqIO;
use Bio::PrimarySeq;  

$dsn = "dbi:mysql:database=simon_oct07;host=ensrv3.sanger.ac.uk";
$db = DBI->connect( $dsn, 'ensadmin' );

open (AGP, shift) or die "couldnt open agp file";

my $file=shift;

while( <AGP> ) {
  split;
  push( @nt, [$_[2], $_[3], $_[6], $_[7],$_[8], $_[4] ] );
  $id = $_[0];
  # global start, global end, clone.version, clone start, clone end, strand
}

my %clones;
foreach my $contig ( @nt ) {
  my ($clonename) = ( $contig->[2] =~ /^([^.]+)/ );
  # print STDERR " Clonename $clonename\n";
  $clones{$clonename} = 1;
}


# fetch all contig and dna information from the db
foreach my $cloneid ( keys %clones ) {
  my $sth = $db->prepare( "select contig.offset, contig.length, dna.sequence from contig, clone, dna where contig.clone = clone.internal_id and clone.id = '$cloneid' and dna.id = contig.dna" );
  $sth->execute;
  while( $arr = $sth->fetchrow_arrayref ) {
    my ( $offset, $length, $seq ) = @$arr;
    push( @dnadata, [ $cloneid, $offset, $length, $seq ] );
  }
}

# now try build in $seq the sequence.
@sortnt = sort { $a->[0] <=> $b->[0] } @nt;
my $seq = "";

#foreach my $dna (@dnadata ) {
#  print join( " ", @{$dna}[0..2] ),"\n";
#}



SEG: foreach my $seg ( @sortnt ) {
  if( $seg->[1]-$seg->[0]-$seg->[4]+$seg->[3] != 0 ) {
    print STDERR join( " ",@$seg ),"\n";
    die "AGP line not reasonable".$seg->[2]."\n";
  }

  my ($clonename) = ( $seg->[2] =~ /^([^.]+)/ );
  foreach my $contig ( @dnadata ) {
    next, unless $clonename eq $contig->[0];
    if( $seg->[3] >= $contig->[1] ) {
      if( $seg->[4] < $contig->[1]+$contig->[2] ) {
        # start end in segment within contig!
	my $start = $seg->[3] - $contig->[1];
	my $len = $seg->[1]-$seg->[0]+1;
	my $subseq = substr( $contig->[3], $start, $len );
	if( $seg->[5] eq "-1" ) {
	  $subseq = reverse( $subseq );
	  $subseq =~ tr/ATCG/TAGC/;
	}
	$seq .= $subseq;
	next SEG;          
      }
    }
  }
  # ups shouldnt be here
  die "Error, shouldnt be here.\n";
}
my $seqout= Bio::SeqIO->new(-file => ">$file" , '-format' => 'Fasta');
my $bioseq=Bio::PrimarySeq->new ( -seq => $seq,
			    -id  => 'ref|NT_000061|Hs2_401 Homo sapiens 2q37.3 sequence /len=35491',
			    -accession => 'NT_000061',
			    -moltype => 'dna'
			    );
$seqout->write_seq($bioseq);


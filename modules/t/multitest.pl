#
# print some information about the example genome
#

use MultiTestDB;


my $multi = MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );

my ($chr, $chr_start, $chr_end);
my ( $contig_id, $contig_name, $id, $name );

my $sth = $db->prepare( "
  select chr.name, min( a.chr_start ), max( a.chr_end )
    from contig c, assembly a, chromosome chr 
   where c.contig_id = a.contig_id
     and chr.chromosome_id = a.chromosome_id
group by a.chromosome_id
" );

$sth->execute();
$sth->bind_columns( \$chr, \$chr_start, \$chr_end );
print "DNA in the assembly\n";
while( $sth->fetch() ) {
  print "   Chr: $chr Start: $chr_start End: $chr_end\n";
}



print "Contig ids and names\n";
$sth = $db->prepare( "select contig_id, name from contig" );
$sth->execute();
$sth->bind_columns( \$contig_id, \$contig_name );
while( $sth->fetch() ) {
  print "   Contig_id: $contig_id Name: $contig_name\n";
}

print "Gene ids and names\n";
$sth = $db->prepare( "select gene_id, stable_id from gene_stable_id" );
$sth->execute();
$sth->bind_columns( \$id, \$name );
while( $sth->fetch() ) {
  print "   dbID: $id Name: $name\n";
}

print "Transcript ids and names.\n";
$sth = $db->prepare( "select transcript_id, stable_id from transcript_stable_id" );
$sth->execute();
$sth->bind_columns( \$id, \$name );
while( $sth->fetch() ) {
  print "   dbID: $id Name: $name\n";
}

print "Translations ids and names\n";
$sth = $db->prepare( "select translation_id, stable_id from translation_stable_id" );
$sth->execute();
$sth->bind_columns( \$id, \$name );
while( $sth->fetch() ) {
  print "   dbID: $id Name: $name\n";
}

print "Analysis ids and names\n";
$sth = $db->prepare( "select * from analysis" );
$sth->execute();
while( my $hr = $sth->fetchrow_hashref() ) {
  print "   ",join( "\n   ", %$hr ),"\n";
  print "---------\n\n";
}



# this script check some of the foreign key relationships in ensembl style 
# database. No commandline args, you have to edit $db= line

use DBI;
use strict; 
use IO::File;

my $db = DBI->connect( "dbi:mysql:host=ecs2d;database=embl_6_29_new", "ensro","" );

# check if contigs are the same in embl and core
my $count = $db->selectrow_array
  ( q{
     select count(*)
       from contig c, homo_sapiens_core_6_29.contig c2 
      where c.id = c2.id 
        and c.offset != c2.offset
        and c.length != c2.length 
      } 
  );

print STDERR "$count differing contigs\n";

# check if new contigs in embl not in core
$count = $db->selectrow_array
  ( q{
       select count(*)
         from contig c
        left join homo_sapiens_core_6_29.contig c2 
          on c.id = c2.id 
         where c2.id is null
  } 
  );
print STDERR "new contigs in embl $count\n";


orphan( $db, "exon", "exon_id", "exon_transcript", "exon_id" );
orphan( $db, "exon_transcript", "exon_id", "exon", "exon_id" );
orphan( $db, "exon", "exon_id", "exon_stable_id", "exon_id" );
orphan( $db, "exon_stable_id", "exon_id", "exon", "exon_id" );


orphan( $db, "transcript", "transcript_id", "exon_transcript", "transcript_id" );
orphan( $db, "exon_transcript", "transcript_id", "transcript", "transcript_id" );
orphan( $db, "transcript", "transcript_id", "transcript_stable_id", "transcript_id" );
orphan( $db, "transcript_stable_id", "transcript_id", "transcript", "transcript_id" );

orphan( $db, "translation", "translation_id", "transcript", "translation_id" );
orphan( $db, "transcript", "translation_id", "translation", "translation_id" );

orphan( $db, "translation", "translation_id", "translation_stable_id", "translation_id" );
orphan( $db, "translation_stable_id", "translation_id", "translation", "translation_id" );


orphan( $db, "gene", "gene_id", "gene_stable_id", "gene_id" );
orphan( $db, "gene_stable_id", "gene_id", "gene", "gene_id" );

orphan( $db, "gene", "gene_id", "transcript", "gene_id" );
orphan( $db, "transcript", "gene_id", "gene", "gene_id" );

orphan( $db, "objectXref", "xrefId", "Xref", "xrefId" );

# very slow, now xref index on obejctXref
# orphan( $db, "Xref", "xrefId", "objectXref", "xrefId" );


exit;

#$count = $db->selectrow_array
#  ( q{
#  
#  } 
#  );
#print STDERR "$count\n";


sub orphan {
  my ( $db, $table1, $col1, $table2, $col2 ) = @_;
  print STDERR "Checking $table1 against $table2\n";

  my $count = $db->selectrow_array
  ( qq{
       select count(*)
         from $table1
        left join $table2
        on $table1.$col1 = $table2.$col2
        where $table2.$col2 is null
      } 
  );

  if( $count > 0 ) {
    print STDERR "$count orphans on $table1, $table2\n";
  } else {
    print STDERR "Ok\n";
  }

}



use DBI;
use strict; 
use IO::File;
use Data::Dumper;

my $db = DBI->connect( "dbi:mysql:host=ecs1c;database=arne_main_29", "ensro","" );
my $count;


# exon_transcript without exon ?
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

orphan( $db, "object_xref", "xref_id", "xref", "xref_id" );

# very slow, now xref index on obejctXref
# orphan( $db, "Xref", "xrefId", "objectXref", "xrefId" );


exit;

$count = $db->selectall_array
  ( q{
  
  } 
  );
print STDERR "$count\n";


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



#
# print some information about the example genome
#
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Test::MultiTestDB;

if( ! @ARGV ) {
  print <<HELP

  This script produces infos about the testgenome MultiTestDB.
  Options are:  -gene 
                -transcript
                -exon 
		-translation 
                -analysis 
                -assembly 
                -contig

HELP
;
  exit;
}

my ( $print_gene, $print_transcript, $print_exon,
     $print_translation, $print_assembly, $print_contig,
     $print_analysis );

GetOptions
  ( 
   'gene' => \$print_gene,
   'transcript' => \$print_transcript,
   'exon' => \$print_exon,
   'translation' => \$print_translation,
   'analysis' => \$print_analysis,
   'assembly' => \$print_assembly,
   'contig' => \$print_contig
  );
	   

my $sth;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );

if( $print_assembly ) {
  my ($chr, $chr_start, $chr_end);
  my ( $contig_id, $contig_name, $id, $name );
  
  $sth = $db->prepare( "
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
    my $slice = $db->get_SliceAdaptor()->fetch_by_chr_start_end
      ( $chr, $chr_start, $chr_end );
    my $tiles = $slice->get_tiling_path();
    for my $tile ( @$tiles ) {
      print "    ";
      print join( " ", "Contig:", $tile->component_Seq->dbID, 
		  "Start:", $tile->component_start, 
		  "End:", $tile->component_end, 
		  "SliceStart:", $tile->assembled_start, 
		  "SliceEnd:", $tile->assembled_end(), 
		  "Ori:", $tile->component_ori() );
      print "\n";
    }    
  }
}

my ( $contig_id, $contig_name, $id, $name );

if( $print_contig ) {

  print "Contig ids and names\n";
  $sth = $db->prepare( "select contig_id, name from contig" );
  $sth->execute();
  $sth->bind_columns( \$contig_id, \$contig_name );
  while( $sth->fetch() ) {
    print "   Contig_id: $contig_id Name: $contig_name\n";
  }
}

if( $print_gene ) {

  print "Gene ids and names\n";
  $sth = $db->prepare( "select gene_id, stable_id from gene_stable_id" );
  $sth->execute();
  $sth->bind_columns( \$id, \$name );
  while( $sth->fetch() ) {
    print "   dbID: $id Name: $name\n";
  }
}


if( $print_transcript ) {
  print "Transcript ids and names.\n";
  $sth = $db->prepare( "select transcript_id, stable_id from transcript_stable_id" );
  $sth->execute();
  $sth->bind_columns( \$id, \$name );
  while( $sth->fetch() ) {
    print "   dbID: $id Name: $name\n";
  }
}

if( $print_translation ) {
  print "Translations ids and names\n";
  $sth = $db->prepare( "select translation_id, stable_id from translation_stable_id" );
  $sth->execute();
  $sth->bind_columns( \$id, \$name );
  while( $sth->fetch() ) {
    print "   dbID: $id Name: $name\n";
  }
}

if( $print_exon ) {
  print "Exon information\n";
  $sth = $db->prepare( "
     select et.transcript_id, e.exon_id, esi.stable_id, 
            e.contig_id, e.contig_start, e.contig_end, e.contig_strand 
       from exon e, exon_stable_id esi, exon_transcript et
      where e.exon_id = esi.exon_id
        and e.exon_id = et.exon_id    
      order by et.transcript_id, et.rank, e.sticky_rank
" );
  $sth->execute();
  my ( $tr_id, $exon_id, $contig_start, $contig_end, $contig_strand, $contig_id,
       $stable_id );

  $sth->bind_columns( \$tr_id, \$exon_id, \$stable_id, \$contig_id, 
		      \$contig_start, \$contig_end, \$contig_strand );
  while( $sth->fetch() ) {
    print "  ";
    print join( " ", "Transcript id $tr_id",
		"exon_id $exon_id", $stable_id,
		"contig $contig_id", 
		"start $contig_start",
		"end $contig_end",
		"strand $contig_strand" );
    print "\n";
  }
}


if( $print_analysis ) {
  print "Analysis ids and names\n";
  $sth = $db->prepare( "select * from analysis" );
  $sth->execute();
  while( my $hr = $sth->fetchrow_hashref() ) {
    print "   ",join( "\n   ", %$hr ),"\n";
    print "---------\n\n";
  }
}




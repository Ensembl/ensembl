# read the mapping out put file
# have the first free id in the variables

# go through new exons.
# mapped, make a exon_stable_id row
# not mapped, generate a new stable_id

# go through genes with same technique

# go through transcript
# if mapped, retrieve the translation_stable_id

use strict;
use DBI;

my $today = "2002-02-20";

my $sdbh = DBI->connect( "DBI:mysql:host=ecs1d;database=homo_sapiens_core_130", "ensro", "");
my $tdbh = DBI->connect( "DBI:mysql:host=ecs1e;database=ens_NCBI_28", "ensro", "" );

my $starttime = scalar( localtime() );

my $sth;
my $resultfile = "map_130_28.txt";
my %stable_id;

for my $table ( "gene_stable_id", "exon_stable_id", 
		"transcript_stable_id", "translation_stable_id" ) {
  $sth = $sdbh->prepare( "select max( stable_id ) from $table" );
  $sth->execute();
  ( $stable_id{$table} ) = $sth->fetchrow_array();
}

# print join( " ", %stable_id ),"\n";
# exit;


exon_stables();
gene_stables();
trans_stables();


sub exon_stables {
  my %exonTranslate;

  open( RES, $resultfile );
  while( <RES> ) {
    /Exon mapped:/ || next;
    my @cols;
    chomp;
    @cols = split( "\t",$_ ); 
    $exonTranslate{$cols[2]} = $cols[3];
  }

  close( RES );
  
  $sth = $tdbh->prepare( "select distinct(exon_id) from exon" );
  $sth->execute;
  open( EXON, ">exon_stable_id.txt" ) or die;

  while( my ( $exonId ) = $sth->fetchrow_array() ) {
    if( ! exists $exonTranslate{$exonId} ) {
      $stable_id{'exon_stable_id'}++;
      print EXON ( join( "\t",($exonId, $stable_id{'exon_stable_id'},"1\t$today\t$today")),"\n" );
    } else {
      print EXON ( join( "\t",($exonId, $exonTranslate{$exonId},"1\t$today\t$today")),"\n" );
    }
  }
  close( EXON );
}

sub gene_stables {
  my %geneTranslate;

  open( RES, $resultfile );
  while( <RES> ) {
    /Gene mapped:/ || next;
    my @cols;
    chomp;
    @cols = split( "\t",$_ ); 
    $geneTranslate{$cols[2]} = $cols[3];
  }

  close( RES );
  
  $sth = $tdbh->prepare( "select gene_id from gene" );
  $sth->execute;
  
  open( GEN, ">gene_stable_id.txt" ) or die;

  while( my ( $geneId ) = $sth->fetchrow_array() ) {
    if( ! exists $geneTranslate{$geneId} ) {
      $stable_id{'gene_stable_id'}++;
      print GEN ( join( "\t",($geneId, $stable_id{'gene_stable_id'},"1\t$today\t$today")),"\n");
    } else {
      print GEN ( join( "\t",($geneId, $geneTranslate{$geneId},"1\t$today\t$today")),"\n" );
    }
  }
  close ( GEN );
}

sub trans_stables {
  my %transcriptHash;
  my %transcriptTranslation;

  open( RES, $resultfile );
  while( <RES> ) {
    /Transcript mapped:/ || next;
    my @cols;
    chomp;
    @cols = split( "\t",$_ ); 
    $transcriptHash{$cols[2]} = $cols[3];
  }

  close( RES );
  
  $sth = $sdbh->prepare( "select tsi.stable_id, tli.stable_id 
                            from transcript_stable_id tsi,
                                 transcript t, translation_stable_id tli
                           where t.transcript_id = tsi.transcript_id 
                             and t.translation_id = tli.translation_id" );
  $sth->execute;
  
  while( my $arrref = $sth->fetchrow_arrayref() ) {
    $transcriptTranslation{$arrref->[0]} = $arrref->[1];
  }

  $sth = $tdbh->prepare( "select transcript_id, translation_id from transcript" );
  $sth->execute();
  
  open( TSC, ">transcript_stable_id.txt" ) or die;
  open( TSL, ">translation_stable_id.txt" ) or die;

  while( my ( $tsc, $tsl ) = $sth->fetchrow_array() ) {
    if( ! exists $transcriptHash{$tsc} ) {
      $stable_id{'transcript_stable_id'}++;
      $stable_id{'translation_stable_id'}++;

      print TSC ( join( "\t",($tsc, $stable_id{'transcript_stable_id'},"1")),"\n");
      print TSL ( join( "\t",($tsl, $stable_id{'translation_stable_id'},"1")),"\n");
    } else {
      my $stable_id = $transcriptHash{$tsc};
      print TSC ( join( "\t",($tsc, $stable_id,"1")),"\n");
      print TSL ( join( "\t",($tsl, $transcriptTranslation{$stable_id},"1")),"\n");
    }
  }
  close ( TCS );
  close ( TSL );
}

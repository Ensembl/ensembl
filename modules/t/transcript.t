use lib 't';
use strict;
use warnings;
use vars qw( $verbose );

BEGIN { $| = 1;  
	use Test;
	plan tests => 17;
}

use MultiTestDB;
use TestUtils qw( debug test_getter_setter );
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Slice;

my $multi = MultiTestDB->new();

$verbose = 0;

ok( $multi );


my $db = $multi->get_DBAdaptor('core' );

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_chr_start_end("20", 30_249_935, 31_254_640 );

my $genes = $slice->get_all_Genes();

my $translates = 1;
my $utr_trans;

debug( "Checking if all Transcripts translate and transform" );

for my $gene ( @$genes ) {
  for my $trans ( @{$gene->get_all_Transcripts()} ) {
    if( $trans->translate()->seq() =~ /\*./ ) {
      $translates = 0;
      debug( $trans->stable_id()." does not translate." );
      last;
    }
    if( $trans->coding_start() != $trans->start() &&
	$trans->coding_end() != $trans->end() ) {
      $utr_trans = $trans->stable_id();
    }
  }

  $gene->transform();

  for my $trans ( @{$gene->get_all_Transcripts()} ) {
    if( $trans->translate()->seq() =~ /\*./ ) {
      $translates = 0;
      debug( $trans->stable_id()." does not translate." );
      last;
    }
  }
}

debug( "utr Transcript is $utr_trans" );

ok( $translates );


my $ta = $db->get_TranscriptAdaptor();

my $tr = $ta->fetch_by_stable_id( "ENST00000217347" );

ok( $tr );

ok( test_getter_setter( $tr, "external_name", "ECSTR-12" ));
ok( test_getter_setter( $tr, "external_db", "EMBL" ));
ok( test_getter_setter( $tr, "dbID", 100000 ));
ok( test_getter_setter( $tr, "type", "NOVEL" ));
ok( $tr->translation->isa( "Bio::EnsEMBL::Translation" ));

debug( "start() == ".$tr->start() );
ok( $tr->start() == 79874 );

debug( "end() == ".$tr->end() );
ok( $tr->end() == 110306 );

debug( "spliced_seq->substr == \"".substr( $tr->spliced_seq(),0, 10 )."\"" );
ok( substr( $tr->spliced_seq(), 0, 10 ) eq "ACGAGACGAA" ); 

debug( "translateable_seq->substr == \"".substr( $tr->translateable_seq(),0,10 )."\"" );
ok( substr( $tr->translateable_seq(),0,10 ) eq "ATGGCAGTGA" );

debug( "coding_start() == ".$tr->coding_start() );
ok( $tr->coding_start() == 85834 );

debug( "coding_end() == ".$tr->coding_end() );
ok( $tr->coding_end() == 108631 );

debug( "pep2genomic: ".($tr->pep2genomic( 10,20))[0]->start());

my @pepcoords = $tr->pep2genomic( 10, 20 );
ok( $pepcoords[0]->start() == 85861 );

debug( "start Exon: ".$tr->start_Exon->stable_id() );
debug( "end Exon: ".$tr->end_Exon->stable_id() );

debug( "five_prime_utr: ".substr( $tr->five_prime_utr()->seq(), -5 , 5 ));
ok( substr( $tr->five_prime_utr()->seq(), -5, 5) eq "CGAAG" ); 

debug( "three_prime_utr: ".substr( $tr->three_prime_utr()->seq(), -5, 5  ));
ok( substr( $tr->three_prime_utr()->seq(), -5, 5  ) eq "TTCAA");


# add_Exon
# get_all_Exons
# length
# get_all_Introns
# flush_Exons
# five_prime_utr
# three_prime_utr

# get_all_translateable_Exons
# translate
# sort
# pep2genomic
# cdna2genomic
# start_Exon
# end_Exon

use lib 't';
use strict;

BEGIN {
    $| = 1;
    use Test;
    plan tests => 37;
}

END {

}


use Bio::EnsEMBL::Utils::Converter;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::RawContig;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;

# Global constants

my $ens_analysis = new Bio::EnsEMBL::Analysis(
    -dbID => 1,
    -logic_name => 'test_fake_analysis'
);

my $ens_contig = new Bio::EnsEMBL::RawContig();

my $seqFeature1 = new Bio::SeqFeature::Generic(
    -start => 100,
    -end => 200,
    -strand => 1,
    -score => 10
);

my $seqFeature2 = new Bio::SeqFeature::Generic(
    -start => 1234,
    -end => 5678,
    -strand => -1,
    -score => 14
);

my $featurePair1 = new Bio::SeqFeature::FeaturePair(
    -feature1 => $seqFeature1,
    -feature2 => $seqFeature2
);

use Bio::EnsEMBL::Utils::Converter;

&test_SeqFeature;
&test_FeaturePair;
&test_hit;
&test_predictionTranscript;
&test_gene;

# Test for SeqFeature, 10 OKs
sub test_SeqFeature{
    # bio to ens
    #
    # testing the analysis conversion firsti
    my $analysis = $ens_analysis;
    eval { require Bio::Pipeline::Analysis;};

    unless($@){
        $analysis = new Bio::Pipeline::Analysis(
            -logic_name => 'analysis to test converter'
        ); 
    }
    my $converter = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::Generic',
        -out => 'Bio::EnsEMBL::SeqFeature',
        -analysis => $analysis,
        -contig => $ens_contig
    );

    my ($ens_seqFeature1, $ens_seqFeature2) = 
        @{$converter->convert([$seqFeature1, $seqFeature2])};

    ok ref($ens_seqFeature1), 'Bio::EnsEMBL::SeqFeature';
    ok ref($ens_seqFeature2), 'Bio::EnsEMBL::SeqFeature';
    
    ok ref($ens_seqFeature1->analysis), 'Bio::EnsEMBL::Analysis';
    
    ok $seqFeature1->start, $ens_seqFeature1->start;
    ok $seqFeature1->end, $ens_seqFeature1->end;
    ok $seqFeature1->strand, $ens_seqFeature1->strand;
    ok $seqFeature1->score, $ens_seqFeature1->score;
    ok $seqFeature1->source_tag, $ens_seqFeature1->source_tag;
    
    # from ens to bio
    $converter = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::EnsEMBL::SeqFeature',
        -out => 'Bio::SeqFeature::Generic',
    );

    my ($bio_seqFeature1, $bio_seqFeature2) =
        @{$converter->convert([$ens_seqFeature1, $ens_seqFeature2])};

    ok ref($bio_seqFeature1), 'Bio::SeqFeature::Generic';
    ok ref($bio_seqFeature2), 'Bio::SeqFeature::Generic';

    
    ok $seqFeature1->start, $bio_seqFeature1->start;
    ok $seqFeature1->end, $bio_seqFeature1->end;
    ok $seqFeature1->strand, $bio_seqFeature1->strand;
    ok $seqFeature1->score, $bio_seqFeature1->score;
    ok $seqFeature1->source_tag, $bio_seqFeature1->source_tag;
    
}


sub test_FeaturePair{
    
    my $converter = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::FeaturePair',
        -out => 'Bio::EnsEMBL::RepeatFeature',
        -analysis => $ens_analysis,
    );
    $converter->contig($ens_contig);
    my ($ens_repeatFeature1) = @{$converter->convert([$featurePair1])};

    ok ref($ens_repeatFeature1), 'Bio::EnsEMBL::RepeatFeature';
    
    ok $ens_repeatFeature1->start, $featurePair1->feature1->start;
    ok $ens_repeatFeature1->end, $featurePair1->feature1->end;
    ok $ens_repeatFeature1->strand, $featurePair1->feature1->strand;
    ok $ens_repeatFeature1->seqname, $featurePair1->feature1->seq_id;
    ok $ens_repeatFeature1->repeat_consensus->name, $featurePair1->feature2->seq_id;
    
    if($featurePair1->feature1->strand == 1){
        ok $ens_repeatFeature1->hstart, $featurePair1->feature2->start;
        ok $ens_repeatFeature1->hend, $featurePair1->feature2->end;
    }else{
        ok $ens_repeatFeature1->hstart, $featurePair1->feature2->end;
        ok $ens_repeatFeature1->hend, $featurePair1->feature2->start;
    }
        
}

# 14 OKs
sub test_hit {
    my $converter = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::Search::HSP::GenericHSP',
        -out => 'Bio::EnsEMBL::BaseAlignFeature',
        -analysis => $ens_analysis,
        -contig => $ens_contig
    );

    use Bio::SearchIO;
    my $searchio = new Bio::SearchIO(
        -format => 'blast',
        -file => 't/converter.blast'
    );

    my @hsps = ();
    while(my $result = $searchio->next_result){
        while(my $hit = $result->next_hit){
            while(my $hsp = $hit->next_hsp){
                push @hsps, $hsp;
            }
        }
    }

    my @align_features = @{$converter->convert(\@hsps)};
    my $align_feature = shift @align_features;
    ok $align_feature->isa('Bio::EnsEMBL::BaseAlignFeature');
    ok $align_feature->start, 1;
    ok $align_feature->end, 504;
    ok $align_feature->strand, 0;
    ok $align_feature->hstart, 21;
    ok $align_feature->hend, 1532;
    ok $align_feature->cigar_string, '504M';

    $align_feature = shift @align_features;
    ok $align_feature->isa('Bio::EnsEMBL::BaseAlignFeature');
    ok $align_feature->start, 64;
    ok $align_feature->end, 324;
    ok $align_feature->strand, 0;
    ok $align_feature->hstart, 182;
    ok $align_feature->hend, 844;
    ok $align_feature->cigar_string, '29M18I22M11I20MD33M4I22M3I25M5I21MI33MD14M';
    
}

sub test_predictionTranscript{
    my $converter = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::Tools::Prediction::Gene',
        -out => 'Bio::EnsEMBL::PreditionTranscript'
    );
    ok 1;
}

sub test_gene{
    my $converter = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::Gene::GeneStructure',
        -out => 'Bio::EnsEMBL::Gene'
    );
    ok 1;
}

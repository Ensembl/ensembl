use lib 't';
use strict;

BEGIN {
    $| = 1;
    use Test;
    plan tests => 22;
}

use Bio::EnsEMBL::Utils::Converter;

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

&test_SeqFeature;
&test_FeaturePair;
&test_hit;

# Test for SeqFeature, 10 OKs
sub test_SeqFeature{
    # bio to ens
    my $converter = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::Generic',
        -out => 'Bio::EnsEMBL::SeqFeature',
        -analysis => $ens_analysis,
        -contig => $ens_contig
    );

    my ($ens_seqFeature1, $ens_seqFeature2) = 
        @{$converter->convert([$seqFeature1, $seqFeature2])};

    ok ref($ens_seqFeature1), 'Bio::EnsEMBL::SeqFeature';
    ok ref($ens_seqFeature2), 'Bio::EnsEMBL::SeqFeature';
    
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
        -contig => $ens_contig
    );

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


sub test_hit {


}


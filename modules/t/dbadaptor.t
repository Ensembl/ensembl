

BEGIN {
    $| = 1;
    $^W = 1;
    print "1..17\n"; 
    use vars qw($loaded);
}

END {
    print "not ok 1\n" unless $loaded;
}

use lib 't';
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use EnsTestDB;

$loaded = 1;

print "ok 1\n";    # 1st test passes.

my( $ens_test_db, $dba );
eval{
    $ens_test_db = EnsTestDB->new;
    $dba = $ens_test_db->get_DBSQL_Obj;
    my $class = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
    die "I was expecting a '$class'"    
        unless ref($dba) eq $class;
};
if ($@) {
    print "not ok 2\n";
    die $@;
} else {
    print "ok 2\n";
}

my $t = 2;

# All these methods contain a "require" statement.
# Check that they all work.
foreach my $dba_get_method (qw{
    get_Feature_Obj
    get_MetaContainer
    gene_Obj
    get_CloneAdaptor
    get_Protfeat_Adaptor
    get_GeneAdaptor
    get_ExonAdaptor
    get_TranscriptAdaptor
    get_FeatureAdaptor
    get_RawContigAdaptor
    get_AnalysisAdaptor
    get_DBEntryAdaptor
    get_StaticGoldenPathAdaptor
    get_KaryotypeBandAdaptor
    get_ChromosomeAdaptor
    })
{
    $t++;
    my( $adaptor );
    eval{
        $adaptor = $dba->$dba_get_method();
    };
    if ($@) {
        warn $@;
        print "not ok $t\n";
    }
    elsif (! $adaptor) {
        warn "No adaptor returned from '$dba_get_method'";
        print "not ok $t\n";
    } else {
        print "ok $t\n";
    }
}

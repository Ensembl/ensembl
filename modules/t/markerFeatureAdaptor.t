use lib 't';
use strict;

BEGIN { $| = 1;  
	use Test ;
	plan tests => 7
}

use MultiTestDB;
use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );

######
# 1  #
######

#test get_MarkerFeatureAdaptor

my $mfa = $db->get_MarkerFeatureAdaptor;

ok($mfa && ref $mfa);

#######
# 2-3 #
#######

#test fetch by dbID

my $mf = $mfa->fetch_by_dbID(1);

ok($mf && ref $mf && $mf->isa('Bio::EnsEMBL::Map::MarkerFeature'));

ok($mf->contig->dbID == 339816 && 
   $mf->analysis->dbID == 10 &&
   $mf->start == 5769 &&
   $mf->end == 5959 &&
   $mf->strand == 0 &&
   $mf->map_weight == 1);


#####
# 4 #
#####


# test fetch all by slice

my $slice = $db->get_SliceAdaptor->fetch_by_chr_start_end('20',
							  '30249935', 
							  '31254640');

my $feats = $mfa->fetch_all_by_Slice($slice);

debug("got " . scalar(@$feats) . " markerfeatures from slice");

ok(scalar(@$feats) == 97);

if($verbose) {
  foreach my $f (@$feats) {
    my $marker = $f->marker;
    print $f->start."-".$f->end ."(". $f->strand.") [".$mf->dbID ."]\n";
    foreach my $ms (@{$marker->get_all_MarkerSynonyms}) {
      print "\t" . $ms->name ."\n";
    }
  }
}

#####
# 5 #
#####

#test fetch_all_by_Marker

my $marker = $db->get_MarkerAdaptor->fetch_by_dbID(2);

$feats = $mfa->fetch_all_by_Marker($marker); 

ok(scalar(@$feats) == 1 && $feats->[0]->start==12671 &&
   $feats->[0]->dbID() == 2);

#######
# 6-7 #
#######

# test store

#hide the contents of the marker_feature table
$multi->hide('core', 'marker_feature');

$marker = $db->get_MarkerAdaptor->fetch_by_dbID(80);
my $contig = $db->get_RawContigAdaptor->fetch_by_dbID(317101);
my $analysis = $db->get_AnalysisAdaptor->fetch_by_dbID(10);

my $marker_feature = Bio::EnsEMBL::Map::MarkerFeature->new;

$marker_feature->contig($contig);
$marker_feature->start(123);
$marker_feature->end(200);
$marker_feature->marker($marker);
$marker_feature->analysis($analysis);

$mfa->store($marker_feature);

ok($marker_feature->dbID && 
   $marker_feature->adaptor == $mfa);


my $sth = $db->prepare('SELECT count(*) from marker_feature');
$sth->execute;
my ($count) = $sth->fetchrow_array;
$sth->finish();

ok($count == 1);


#restore the marker feature table to its original state
$multi->restore('core', 'marker_feature');


#!/usr/local/bin/perl

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Getopt::Long;

my $host      = 'ecs1f';
my $dbname    = 'alistair_mouse_si_Nov01';
my $newhost   = 'ecs1f';
#my $newhost  = 'ecs1e';
my $newdbname = 'mouse_Sanger_Nov01_denormalised';
#my $newdbname = 'mouse_tmp_test';
my $newpath   = 'CHR';
my $dbuser    = 'ensadmin';
my $pass      = 'ensembl';
my $path      = 'SI_Nov01';
my $feature   = 0;
my $start     = 0;
my $end       = 0;
my $chrname;

$| = 1;

&GetOptions( 'host:s'    => \$host,
             'dbuser:s'  => \$dbuser,
             'dbname:s'  => \$dbname,
	     'feature'   => \$feature,
             'path=s'    => \$path,
             'pass:s'    => \$pass,
	     'start:n'   => \$start,
	     'end:n'   => \$end,
	     'chrname:s'   => \$chrname,
             );

print STDERR "Time in dumpgff " . time . "\n";

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host             => $host,
					    -user             => $dbuser,
					    -dbname           => $dbname,
                                            -pass             => $pass,
					    -perlonlyfeatures => 0);

my $newdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host             => $newhost,
					       -user             => $dbuser,
					       -dbname           => $newdbname,
					       -pass             => $pass,
					       -perlonlyfeatures => 0);

print STDERR "Time after db  " . time . "\n";

print STDERR "Connected to database\n";

$db->static_golden_path_type($path);
$newdb->static_golden_path_type($newpath);

my $stadp = $db->get_StaticGoldenPathAdaptor();

my $from_vc = $db->get_StaticGoldenPathAdaptor->fetch_VirtualContig_by_chr_start_end($chrname,
										     $start,
										     $end
										    );

my $to_vc   = $newdb->get_StaticGoldenPathAdaptor->fetch_VirtualContig_by_chr_start_end($chrname,
											$start,
											$end
										       );	

#print STDERR "old vc seq:" . $from_vc->primary_seq->seq . "\n";

my @simfeat     = $from_vc->get_all_SimilarityFeatures;

print STDERR "Got " . scalar(@simfeat) . "similarity features\n";

# now for fun with genscans
# Virtual::Contig->get_all_Prediction features cannot map from RC to VC - "unable to map" errors
# Virtual->StaticContig->get_all_PredictionFeatures returns genscans as Seqfeatures not FeaturePairs so is useless here.
# so ... we'll have to filter them out of the similarity featuers and do wicked things to them to avoid off by 1
# in fact, we ought to do the off by 1 correction for all the similarity features
my @genscans;
foreach my $f(@simfeat){
  if($f->analysis->db eq 'HumanIso.smat'){
    $f->start($f->start + 1);
    $f->end($f->end + 1);
    push(@genscans,$f);
  }
}
print STDERR "Got " . scalar(@genscans) . " genscan features\n";

# get simple features here. We do it by type because otherwise we'll get all the 
# similarity features in too and will have to filter them out. waste of time.
my @tRNA = $from_vc->get_all_SimpleFeatures_by_feature_type("tRNA");
print STDERR "Got " . scalar(@tRNA) . " tRNA features\n";

my @cpg = $from_vc->get_all_SimpleFeatures_by_feature_type("cpg_island");
print STDERR "Got " . scalar(@cpg) . " cpg features\n";

# fake up simple features as FeaturePairs and add them into @allfeat. 
# YES! I know this is mad but it's the easiest way till we switch schemas.
my @allfeat;
push(@allfeat, @simfeat);
push(@allfeat, @genscans);
my @simplefeat;
push(@simplefeat, @tRNA);
push(@simplefeat, @cpg);

foreach my $feat(@simplefeat){
  my $fp = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
  $fp->set_all_fields(
		      $feat->start,
		      $feat->end, 
		      $feat->strand,
		      $feat->score,
		      $feat->primary_tag,
		      $feat->source_tag,
		      'wibble',
		      -1,
		      -1,
		      1,
		      0,
		      $feat->primary_tag,
		      $feat->source_tag,
		      '__NONE__',
		      'NULL',
		      0,
		      0,
		      0
		     );
  $fp->analysis($feat->analysis);
  push (@allfeat, $fp);
}

# features are in chromosomal coords; write them to newdb
# need to convert to raw contig coords. Probably faster to do a straight insert than to use FeatureAdaptor.


# may want to introduce some filtering here - not putting through v.low score blast, eg?
my $converted;
FEAT:
foreach my $f (@allfeat) {
  # convert to denormalised contig coords
  my ($scontig, $fstart, $sstrand) = $to_vc->_vmap->raw_contig_position($f->start, $f->strand);
  my ($econtig, $fend, $estrand) = $to_vc->_vmap->raw_contig_position($f->end, $f->strand);
  
  if($econtig->id != $scontig->id){
    print STDERR "sticky feature - can't map it\n";
    next;
#    return undef;
  }

  if($estrand != $sstrand){
    print STDERR "strand mismatch - can;t map it\n";
    next;
#    return undef;
  }


# print to a tab delimited file that can later be loaded up
  my $query = "insert into feature values(\\N," . $scontig->internal_id . ", $fstart, $fend, " . $f->score . ", $sstrand, " . $f->analysis->dbID . ", '" . $f->analysis->program . "', " . $f->hstart . ", " . $f->hend . ", '" . $f->hseqname . "', '" . $f->p_value . "', " . $f->percent_id . ", " . $f->phase . ", " . $f->end_phase . ")";
  
  print "$query\n";  
  
  # or you could uncomment these two lines to allow direct insertion
  #my $sth =  $newdb->prepare($query);
  
  #my $res = $sth->execute;

}
    


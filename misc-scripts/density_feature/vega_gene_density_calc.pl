#
# script to calculate the gene density features on a database
# should work on any species database
#

#
# It will only run on databases with genes ...
# boundary condition: on average there should be 2 genes per block
#

# this is a modified version of gene_density.calc.pl written for vega by st3 - is a bit slow since when looking for each type of gene it actually retrieves all genes but then only records those of the correct type.


use strict;

use lib '../../modules/','../../../bioperl-live';
use Data::Dumper;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname  );

my ( $block_count, $genome_size, $block_size );

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname
	  );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);

my $sth = $db->prepare( "select count(*) from gene" );
$sth->execute();

my ( $gene_count )  = $sth->fetchrow_array();

if( ! $gene_count ) {
  print STDERR "No gene density for $dbname.\n";
  exit();
} else {
  $block_count = $gene_count >> 1;
}

#
# Could be database without seq_regions
#  Then have to try and attach core db
#
$sth = $db->prepare( "select count(*)  from seq_region" );
$sth->execute();
my ( $seq_region_count ) = $sth->fetchrow_array();
if( ! $seq_region_count ) {
  #
  # for the time being only do core dbs
  # no dbs with no seq regions
  #
  print STDERR "No gene density for $dbname, no seq_regions.\n";
  exit();

  if( ($dbname =~ /_estgene_/ ) || ( $dbname =~ /_vega_/ )) {
    my $dna_db_name = $dbname;
    $dna_db_name =~ s/(_estgene_|_vega_)/_core_/;
    my $dna_db =  new Bio::EnsEMBL::DBSQL::DBAdaptor
      (-host => $host,
       -user => $user,
       -port => $port,
       -pass => $pass,
       -dbname => $dna_db_name );
    print STDERR "Attaching $dna_db_name to $dbname.\n";
    $db->dnadb( $dna_db );
  } else {
    print STDERR "No gene density for $dbname, no seq_regions.\n";
    exit();
  }
}

#
# Get the adaptors needed;
#

my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

#
# block size estimation
#

my $top_slices = $slice_adaptor->fetch_all('toplevel');
for my $slice ( @$top_slices ) {
  $genome_size += $slice->length();
}

$block_size = int( $genome_size / $block_count );

#which types of geens are we interested in?
my $analysis_names = {
		      'pseudoGeneDensity'=>[
					    'Processed_pseudogene',
					    'Pseudogene',
					    'Unprocessed_pseudogene'],
		      'knownGeneDensity'=>'Known',
		      'novelCDSDensity'=>'Novel_CDS',
		      'novelTransDensity'=>'Novel_Transcript',
		      'putativeTransDensity'=>'Putative',
		      'predictedTransDensity'=>'Predicted_Gene',
		      'predictedIgPseudoDensity'=>'Ig_Pseudogene_Segment',
		      'IgSegDensity'=>'Ig_Segment'
		     };
my $total_successes;
my $known_successess;
my $novel_cds_successes;
my $novel_trans_successes;
my $putative_successes;
my $pseudo_successes;
my $predicted_trans_successes;
my $predictedIGPseudo_successes;
my $IgSegment_successes;
my $failures;

# Now the actual feature calculation loop, looking for each gene type in turn

foreach my $anal_name (keys %$analysis_names) {

    my %gene_names;

#create analysis object for each gene type
    my $analysis = new Bio::EnsEMBL::Analysis (-program     => "vega_gene_density_calc.pl",
					    -database    => "ensembl",
					    -gff_source  => "vega_gene_density_calc.pl",
					    -gff_feature => "density",
					    -logic_name  => $anal_name);
    $aa->store( $analysis );

    print "FETCHING DATA ON GENES OF TYPE $analysis_names->{$anal_name}:\n"; # I know this doesn't display correctly for pseudogenes, tough!

    $analysis = $aa->fetch_by_logic_name($anal_name);

  #
  # Create new density type for each gene type
    my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					  -block_size => $block_size,
					  -value_type => 'sum');
    $dta->store($dt);

    my ( $current_start, $current_end );
    
    foreach my $slice (@$top_slices){
	$current_start = 1;
	my @density_features=();
	print "Gene densities for ".$slice->seq_region_name().
	" with block size $block_size\n";
	
	while($current_start <= $slice->end()) {
	    $current_end = $current_start+$block_size-1;
	    if( $current_end > $slice->end() ) {
		$current_end = $slice->end();
	    }
	    
	    my $sub_slice = $slice->sub_Slice( $current_start, $current_end );
	    my $count =0;
	    
      #
      # retrieve info for all genes on the subslice, but count only those of a particular type
      #

	    foreach my $gene (@{$sub_slice->get_all_Genes()}){
#only count each gene once for each gene type loop
		if (exists $gene_names{$gene->stable_id}) {
		    next;
		}
		else {
		    $gene_names{$gene->stable_id} = 1;
		}
#pseudogenes first...
		my $classification = $gene->type;
		if ($anal_name eq 'pseudoGeneDensity'){
		    my $check = 0;
		    foreach my $pseudo_type (@{$analysis_names->{$anal_name}}) {
			if (($classification eq $pseudo_type) & ($check<1)) {
			    $count++;
			    $check = 1;
			    $total_successes++;
			    $pseudo_successes++;
			}
			else {
			    $failures++; 
			}
		    }
		}
#then the other types 
		elsif ($classification eq $analysis_names->{$anal_name}) {
		    if ($anal_name eq 'knownGeneDensity') {
			$count++;
			$total_successes++;
			$known_successess++;
		    } elsif ($anal_name eq 'novelCDSDensity') {
			$count++;
			$total_successes++;
			$novel_cds_successes++;
		    } elsif ($anal_name eq 'novelTransDensity') {
			$count++;
			$total_successes++;
			$novel_trans_successes++;
		    } elsif ($anal_name eq 'putativeTransDensity') {
			$count++;
			$total_successes++;
			$putative_successes++;
		    } elsif ($anal_name eq  'predictedTransDensity') {
			$count++;
			$total_successes++;
			$predicted_trans_successes++;
		    } elsif ($anal_name eq 'predictedIgPseudoDensity') {
			$count++;
			$total_successes++;
			$predictedIGPseudo_successes++;
		    } elsif ($anal_name eq 'IgSegDensity') {
			$count++;
			$total_successes++;
			$IgSegment_successes++;
		    }
		    else {
			$failures++;
		    }
		}
	    }
	    
	    push @density_features, Bio::EnsEMBL::DensityFeature->new
	    (-seq_region    => $slice,
	     -start         => $current_start,
	     -end           => $current_end,
	     -density_type  => $dt,
	     -density_value => $count);
	    
	    $current_start = $current_end + 1;
	    #      print STDERR ".";
	}
	$dfa->store(@density_features);
	print "Created ", scalar @density_features, " gene density features.\n";
	#    print_features(\@density_features);
    }
}
print "number of genes in db =  $gene_count\n";
print "total number of succesfull identifications = $total_successes\n";
print "number of known genes = $known_successess\n";
print "number of novel cds = $novel_cds_successes\n";
print "number of novel trans = $novel_trans_successes\n";
print "number of putative trans = $putative_successes\n";
print "number of pseudogenes = $pseudo_successes\n";
print "number of predicted transgenes = $predicted_trans_successes\n";
print "number of IG Pseudogene Segments = $predictedIGPseudo_successes\n";
print "number of IG Segments = 	$IgSegment_successes\n";
# this gives totally the wrong number of failures and is not to be relied upon at all!
#if ($failures > 0) {
#    print "warning, you have some unmatched genes: a lot less than $failures\n";
#}

#
# helper to draw an ascii representation of the density features
#
sub print_features {
  my $features = shift;

  return if(!@$features);

  my $sum = 0;
  my $length = 0;
#  my $type = $features->[0]->{'density_type'}->value_type();

  print("\n");
  my $max=0;
  foreach my $f (@$features) {
    if($f->density_value() > $max){
      $max=$f->density_value();
    }
  }
  if( !$max ) { $max = 1 };

  foreach my $f (@$features) {
    my $i=1;
    for(; $i< ($f->density_value()/$max)*40; $i++){
      print "*";
    }
    for(my $j=$i;$j<40;$j++){
      print " ";
    }
    print "  ".$f->density_value()."\t".$f->start()."\n";
  }
}




  



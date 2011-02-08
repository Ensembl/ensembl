#!/usr/local/ensembl/bin/perl -w

#
# script to calculate the gene density features on a database
# should work on any species database
#

#
# It will only run on databases with genes ...
#

# I think the right thing here is to generate densities on the longest
# 125 toplevel slices... The website will be happy with about 150 bins I
# think.


use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Getopt::Long;

my $bin_count  = 150;
my $max_slices = 100;

my ( $host, $user, $pass, $port, $dbname, $pattern  );

$port = 3306 ; 

my ( $block_count, $genome_size, $block_size );

GetOptions(
  "host=s", \$host,
  "user=s", \$user,
  "pass=s", \$pass,
  "port=i", \$port,
  "dbname=s", \$dbname,
  "pattern=s", \$pattern,
);

unless ($host || $user || $pass || $dbname || $pattern) {
  print "\n\nusage : perl gene_density.pl -host <HOST> -user <USER> -pass <PASS> -port <3306> -dbname <DATABASENAME>|-pattern <PATTERN> \n\n" ;
  exit(0) ;
}
my @dbnames;
if (! $dbname) {
  my $dsn = sprintf( 'dbi:mysql:host=%s;port=%d', $host, $port );
  my $dbh = DBI->connect( $dsn, $user, $pass );
  @dbnames =
    map { $_->[0] } @{ $dbh->selectall_arrayref('SHOW DATABASES') };
}
else {
  @dbnames = ( $dbname );
}
  
foreach my $dbname (@dbnames) {
  if ( $pattern && ($dbname !~ /$pattern/) ) { next }

  printf( "Connecting to '%s'\n", $dbname );
  
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $host,
    -user => $user,
    -port => $port,
    -pass => $pass,
    -dbname => $dbname);
  
  my $sth = $db->dbc()->prepare( "select count(*) from gene" );
  $sth->execute();

  my ( $gene_count )  = $sth->fetchrow_array();

  if( ! $gene_count ) {
    print STDERR "No gene density for $dbname.\n";
    exit();
  }

#
# Could be database without seq_regions
#  Then have to try and attach core db
#
  $sth = $db->dbc()->prepare( "select count(*)  from seq_region" );
  $sth->execute();
  my ( $seq_region_count ) = $sth->fetchrow_array();
  if( ! $seq_region_count ) {
  #
  # for the time being only do core dbs
  # no dbs with no seq regions
  #
    print STDERR "No gene density for $dbname, no seq_regions.\n";
    exit();

    if( ($dbname =~ /_est_/ ) || ( $dbname =~ /_vega_/ ) || ( $dbname =~ /_cdna_/ ) ) {
      my $dna_db_name = $dbname;
      $dna_db_name =~ s/(_estgene_|_vega_|_cdna_)/_core_/;
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


  print "Deleting old knownGeneDensity and geneDensity features\n";
  $sth = $db->dbc->prepare("DELETE df, dt, a, ad FROM density_feature df, density_type dt, analysis a, analysis_description ad WHERE ad.analysis_id = a.analysis_id AND a.analysis_id=dt.analysis_id AND dt.density_type_id=df.density_type_id AND a.logic_name IN ('knowngenedensity', 'genedensity')");
  $sth->execute();
  
  $sth = $db->dbc->prepare("DELETE df, dt, a FROM density_feature df, density_type dt, analysis a WHERE a.analysis_id=dt.analysis_id AND dt.density_type_id=df.density_type_id AND a.logic_name IN ('knowngenedensity', 'genedensity')");
  $sth->execute();
  

# $sth = $db->dbc()->prepare(
#   qq(
#   DELETE ad
#   FROM analysis_description ad
#   WHERE ad.display_label IN ('knownGeneDensity', 'geneDensity')) );
# $sth->execute();

  my $dfa = $db->get_DensityFeatureAdaptor();
  my $dta = $db->get_DensityTypeAdaptor();
  my $aa  = $db->get_AnalysisAdaptor();
  my $slice_adaptor = $db->get_SliceAdaptor();
  
  
  # Sort slices by coordinate system rank, then by length.
  my @sorted_slices =
    sort( {      $a->coord_system()->rank() <=> $b->coord_system()->rank()
		   || $b->seq_region_length() <=> $a->seq_region_length()
		 } @{ $slice_adaptor->fetch_all('toplevel') } );
  
  my $analysis =
    new Bio::EnsEMBL::Analysis(
      -program     => "gene_density_calc.pl",
      -database    => "ensembl",
      -gff_source  => "gene_density_calc.pl",
      -gff_feature => "density",
      -logic_name  => "knowngenedensity",
      -description => 'Known gene density features in a database ' . 'as calculated by gene_density_calc.pl',
      -display_label => 'Genes (density)',
      -displayable   => 1 );
  
  $aa->store($analysis);
  $aa->update($analysis);

  $analysis = new Bio::EnsEMBL::Analysis (
    -program     => "gene_density_calc.pl",
    -database    => "ensembl",
    -gff_source  => "gene_density_calc.pl",
    -gff_feature => "density",
    -logic_name  => "genedensity",
    -description => 'Gene density features in a database ' . 'as calculated by gene_density_calc.pl',
    -display_label => 'Genes (density)',
    -displayable   => 1 );
  
  $aa->store( $analysis );
  $aa->update($analysis);

#
# Now the actual feature calculation loop
#


  foreach my $known (1, 0) {
  #
  # Create new analysis object for density calculation.
  #

    if($known) {
      $analysis = $aa->fetch_by_logic_name('knowngenedensity');
    } else {
      $analysis = $aa->fetch_by_logic_name('genedensity');
    }

  #
  # Create new density type.
  #

    my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					  -region_features => $bin_count,
					  -value_type => 'sum');

    $dta->store($dt);

    my $slice_count = 0;
    my ( $current, $current_start, $current_end  );
  
    foreach my $slice (@sorted_slices){

      $block_size = $slice->length() / $bin_count;

      my @density_features;#sf7

      print "Gene densities for ".$slice->seq_region_name().
	" with block size $block_size\n";
      $current_end = 0;
      $current = 0;

      while($current_end < $slice->length) {
	$current += $block_size;
	$current_start = $current_end+1;
	$current_end = int( $current + 1 );
	
	if( $current_end < $current_start ) { 
	  $current_end = $current_start;
	}
	
	if( $current_end > $slice->end ) {
	  $current_end = $slice->end;
	}
	

	my $sub_slice = $slice->sub_Slice( $current_start , $current_end );
	
	my $count =0;

      #
      # Store info for genes (ignore pseudo genes)
      #

	foreach my $gene (@{$sub_slice->get_all_Genes()}){
	  if($gene->biotype() !~ /pseudogene/i and $gene->start >=1 ) {
	    $count++ if(!$known || $gene->is_known());
	  }
	}

	push @density_features, Bio::EnsEMBL::DensityFeature->new
	  (-seq_region    => $slice,
	   -start         => $current_start,
	   -end           => $current_end,
	   -density_type  => $dt,
	   -density_value => $count);
      }
      
      $dfa->store(@density_features);
      print "Created ", scalar @density_features, " gene density features.\n";

      last if ( $slice_count++ > $max_slices );
    }
  }
  print "Finished with $dbname";
}


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




  



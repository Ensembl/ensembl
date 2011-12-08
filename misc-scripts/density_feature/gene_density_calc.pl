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
use Bio::EnsEMBL::Utils::ConversionSupport;

my $bin_count  = 150;
my $max_slices = 100;

my ( $host, $user, $pass, $port, $dbname, $pattern  );

my ( $mhost, $mport ) = ( 'ens-staging1', '3306' );
my ( $muser, $mpass ) = ( 'ensro',        undef );
my $mdbname = 'ensembl_production';

$port = 3306 ; 

my ( $block_count, $genome_size, $block_size );

GetOptions(
  "host|h=s", \$host,
  "user|u=s", \$user,
  "pass|p=s", \$pass,
  "port=i",   \$port,
  "dbname|d=s", \$dbname,
  "pattern=s", \$pattern,
            "mhost=s", \$mhost,
            "mport=i", \$mport,
            "muser=s", \$muser,
  "help" ,               \&usage
);

usage() if (!$host || !$user || !$pass || (!$dbname && !$pattern));

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

  print STDOUT "Connecting to $dbname\n";
  
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
      print STDOUT "Attaching $dna_db_name to $dbname.\n";
      $db->dnadb( $dna_db );
    } else {
      print STDERR "No gene density for $dbname, no seq_regions.\n";
      exit();
    }
  }

#
# Get the adaptors needed;
#


  print STDOUT "Deleting old knownGeneDensity and geneDensity features\n";
  
  $sth = $db->dbc->prepare("DELETE df, dt FROM density_feature df, density_type dt, analysis a WHERE a.analysis_id=dt.analysis_id AND dt.density_type_id=df.density_type_id AND a.logic_name IN ('knowngenedensity', 'genedensity')");
  $sth->execute();
  
  my $dfa = $db->get_DensityFeatureAdaptor();
  my $dta = $db->get_DensityTypeAdaptor();
  my $aa  = $db->get_AnalysisAdaptor();
  my $slice_adaptor = $db->get_SliceAdaptor();
  
  my $support = 'Bio::EnsEMBL::Utils::ConversionSupport';
  my $analysis1 = $aa->fetch_by_logic_name('knowngenedensity');
  my $analysis2 = $aa->fetch_by_logic_name('genedensity');

  # Master database location:
  my ( $muser, $mpass ) = ( 'ensro',        undef );
  my $mdbname = 'ensembl_production';
  my $prod_dsn;
  my $prod_dbh;

  if ( !defined($analysis1) || !defined($analysis2) ) {

    #get analyses descriptions from the master database

     $prod_dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $mhost, $mport, $mdbname );
     $prod_dbh = DBI->connect( $prod_dsn, $muser, $mpass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );
  } 

  if ( !defined($analysis1) ) {

   my ($display_label,$description) = $prod_dbh->selectrow_array("select distinct display_label, description from analysis_description where is_current = 1 and logic_name = 'knowngenedensity'");

    $analysis1 = new Bio::EnsEMBL::Analysis(
    -program     => "gene_density_calc.pl",
    -database    => "ensembl",
    -gff_source  => "gene_density_calc.pl",
    -gff_feature => "density",
    -logic_name  => "knowngenedensity",
    -description => $description,
    -display_label => $display_label,
    -displayable   => 1 );

    $aa->store($analysis1);

  } else {
      $analysis1->created($support->date());
      $aa->update($analysis1);
  }
 

  if ( !defined($analysis2) ) {

  my ($display_label,$description) = $prod_dbh->selectrow_array("select distinct display_label, description from analysis_description where is_current = 1 and logic_name = 'genedensity'");

      $analysis2 = new Bio::EnsEMBL::Analysis(
    -program     => "gene_density_calc.pl",
    -database    => "ensembl",
    -gff_source  => "gene_density_calc.pl",
    -gff_feature => "density",
    -logic_name  => "genedensity",
    -description => $description,
    -display_label => $display_label,
    -displayable   => 1 );

      $aa->store($analysis2);  

  } else {
      $analysis2->created($support->date());
      $aa->update($analysis2);
  }

  if ( defined($prod_dbh) ) {
      $prod_dbh->disconnect;
  }

#
# Now the actual feature calculation loop
#
  my $total_count = 0;

  foreach my $known (1, 0) {
  #
  # Create new analysis object for density calculation.
  #
    my $analysis;

    if($known) {
      $analysis = $aa->fetch_by_logic_name('knowngenedensity');
    } else {
      $analysis = $aa->fetch_by_logic_name('genedensity');
    }
    
  # Sort slices by coordinate system rank, then by length.
    my @sorted_slices =
	sort( {      $a->coord_system()->rank() <=> $b->coord_system()->rank()
		   || $b->seq_region_length() <=> $a->seq_region_length()
		 } @{ $slice_adaptor->fetch_all('toplevel') } );


  #
  # Create new density type.
  #

    my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
					  -region_features => $bin_count,
					  -value_type => 'sum');

    $dta->store($dt);

    my $slice_count = 0;
    my ( $current, $current_start, $current_end  );
  
    while ( my $slice = shift @sorted_slices){

      $block_size = $slice->length() / $bin_count;

      my @density_features;#sf7

      print STDOUT "Gene densities for ".$slice->seq_region_name().
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

	$total_count++;
      }
      
      $dfa->store(@density_features);

      last if ( $slice_count++ > $max_slices );
    }

  }
  print STDOUT "Created $total_count gene density features.\n";
  print STDOUT "Finished with $dbname";
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

sub usage {
  my $indent = ' ' x length($0);
  print <<EOF; exit(0);


What does it do?

Populates known gene density features in a database as well as gene density
features from genes of any status.

First it needs gene and seq_region tables to be populated. It then
deletes all knownGeneDensity and geneDensity entries from the
analysis, density_type and density_feature tables. All toplevel
slices are fetched and sorted from longest to shortest. Each slice
is divided into 150 bins (blocks). For each of the blocks or
sub-slices, the number of genes is counted to give geneDensity and
the number of genes with status KNOWN status is counted to give
knownGeneDensity. All biotypes except pseudogene are counted.

Input data: genes, top level seq regions, xrefs 
Output tables: analysis (logic_name: knownGeneDensity and geneDensity), 
               analysis description, density_type, density_feature


When to run it in the release cycle?

It can be run after compara have handed over homologies and core have 
done xref projections.


Which databases to run it on?

It needs to be run across all core databases for every release.


How long does it take?

It takes about 10 mins to run for a database in normal queue,


Usage:

  $0 -h host [-port port] -u user -p password \\
  $indent -d database | -pattern pattern \\
  $indent [-help]  \\


  -h|host              Database host to connect to

  -port                Database port to connect to (default 3306)

  -u|user              Database username

  -p|pass              Password for user

  -d|dbname            Database name


  -mhost              ensembl_production database host to connect to

  -mport              ensembl_production database port to connect to

  -muser              ensembl_production database username


  -pattern             Database name regexp

  -help                This message


EOF

}


  



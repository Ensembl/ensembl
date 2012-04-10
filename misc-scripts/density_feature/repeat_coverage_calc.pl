#!/usr/local/ensembl/bin/perl -w
#
# Calculate the repeat coverage for given database.
# condition: 1k blocks to show contigview displays
#  4000 blocks for a whole genome
#
# checks wether database contains repeats before doing anything

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DensityTypeAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Utils::ConversionSupport;

use POSIX;

use strict;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname  );

# Master database location:
my ( $mhost, $mport ) = ( 'ens-staging1', '3306' );
my ( $muser, $mpass ) = ( 'ensro',        undef );
my $mdbname = 'ensembl_production';


GetOptions( "host|h=s", \$host,
	    "user|u=s", \$user,
	    "pass|p=s", \$pass,
	    "port=i",   \$port,
	    "dbname|d=s", \$dbname,
            "mhost=s", \$mhost,
            "mport=i", \$mport,
            "muser=s", \$muser,
	    "help" ,               \&usage
	  );

usage() if (!$host || !$user || !$pass || !$dbname);

my $chunksize = 1_000_000;
my $small_blocksize = 1_000;
my $bin_count = 150;
my $max_top_slice = 100;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);



my $sth = $db->dbc()->prepare( "select count(*) from repeat_feature" );
$sth->execute();

my ( $repeat_count )  = $sth->fetchrow_array();

if( ! $repeat_count ) {
  print STDERR "No repeat density for $dbname.\n";
  exit();
}

#
# Get the adaptors needed;
#

#
# Clean up old features first. Also remove analysis and density type entry as these are recreated.
#

print STDOUT "Deleting old PercentageRepeat features\n";
$sth = $db->dbc->prepare("DELETE df, dt FROM density_feature df, density_type dt, analysis a WHERE a.analysis_id=dt.analysis_id AND dt.density_type_id=df.density_type_id AND a.logic_name='percentagerepeat'");
$sth->execute();


my $slice_adaptor = $db->get_SliceAdaptor();
my $dfa = $db->get_DensityFeatureAdaptor();
my $dta = $db->get_DensityTypeAdaptor();
my $aa  = $db->get_AnalysisAdaptor();


my $analysis = $aa->fetch_by_logic_name('percentagerepeat');


if ( !defined($analysis) ) {


   my $prod_dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $mhost, $mport, $mdbname );
   my $prod_dbh = DBI->connect( $prod_dsn, $muser, $mpass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

   my ($display_label,$description) = $prod_dbh->selectrow_array("select distinct display_label, description from analysis_description where is_current = 1 and logic_name = 'percentagerepeat'");

   $prod_dbh->disconnect;

   $analysis = new Bio::EnsEMBL::Analysis(
              -program     => "repeat_coverage_calc.pl",
              -database    => "ensembl",
              -gff_source  => "repeat_coverage_calc.pl",
              -gff_feature => "density",
              -logic_name  => "percentagerepeat", 
              -description => $description,
              -display_label => $display_label,
              -displayable   => 1 );

    $aa->store($analysis);
} else {

    my $support = 'Bio::EnsEMBL::Utils::ConversionSupport';
    $analysis->created($support->date());
    $aa->update($analysis);

}

my $slices = $slice_adaptor->fetch_all( "toplevel" );
my @sorted_slices = sort { $b->seq_region_length() <=> $a->seq_region_length() } @$slices;

my $small_density_type = Bio::EnsEMBL::DensityType->new
  (-analysis   => $analysis,
   -block_size => $small_blocksize,
   -value_type => 'ratio');

my $variable_density_type = Bio::EnsEMBL::DensityType->new
  (-analysis   => $analysis,
   -region_features => $bin_count,
   -value_type => 'ratio');

$dta->store($small_density_type);
$dta->store($variable_density_type);


my $slice_count = 0;


while ( my $slice = shift @sorted_slices ) {
  
  #
  # do it for small and large blocks
  #
  print STDOUT ("Working on seq_region ".$slice->seq_region_name()." length ".$slice->seq_region_length());

  my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();
  my $chunk_end = 0;
  my $variable_end = 0;
  my $small_end = 0;
  my $small_start = 0;
  my $repeat_size;
  my $variable_start = 0;
  my $variable_blocksize = POSIX::ceil( $slice->seq_region_length() / 
					$variable_density_type->region_features());
  $slice_count++;

  while( $chunk_end < $slice->length() ) {
    my $chunk_start = $chunk_end+1;
    $chunk_end += $chunksize;
    $chunk_end = $slice->length() if $chunk_end > $slice->length();

    register( $rr, $slice, $chunk_start, $chunk_end );

    my @dfs = ();

    if( $slice_count < $max_top_slice ) {
      while ( $variable_end+$variable_blocksize <= $chunk_end ) {
	# here we can do the variable sized repeat densities
	$variable_start = $variable_end+1;
	$variable_end += $variable_blocksize;

	$repeat_size = $rr->overlap_size( "1", $variable_start, $variable_end );
	my $percentage_repeat = $repeat_size / $variable_blocksize * 100;

	push( @dfs, Bio::EnsEMBL::DensityFeature->new
	      (-seq_region    => $slice,
	       -start         => $variable_start,
	       -end           => $variable_end,
	       -density_type  => $variable_density_type,
	       -density_value => $percentage_repeat));

      }
    }

    while ( $small_end + $small_blocksize <= $chunk_end ) {
      # here we can do the small sized density features
      $small_start = $small_end+1;
      $small_end += $small_blocksize;

      $repeat_size = $rr->overlap_size( "1", $small_start, $small_end );
      my $percentage_repeat = $repeat_size / $small_blocksize * 100;

      push( @dfs, Bio::EnsEMBL::DensityFeature->new
        (-seq_region    => $slice,
         -start         => $small_start,
         -end           => $small_end,
         -density_type  => $small_density_type,
         -density_value => $percentage_repeat));
    }

    if (@dfs) {
      $dfa->store(@dfs);
    } else {
      warning("No repeat density calculated for ".$slice->name." (chunk start $chunk_start, chunk end $chunk_end).");
    }

    my $used_lower_limit = $small_start < $variable_start ? $small_start : $variable_start;

    # here some rr cleanup
    $rr->check_and_register( "1", 0, $used_lower_limit );
  }

  # missing the last bits
  if( $small_end < $slice->length() ) {
    $small_start = $small_end+1;
    $small_end = $slice->length();

    $repeat_size = $rr->overlap_size( "1", $small_start, $small_end );
    my $percentage_repeat = $repeat_size / ($small_end - $small_start + 1 ) * 100;

    $dfa->store( Bio::EnsEMBL::DensityFeature->new
		 (-seq_region    => $slice,
		  -start         => $small_start,
		  -end           => $small_end,
		  -density_type  => $small_density_type,
		  -density_value => $percentage_repeat));
  }

  if( $variable_end < $slice->length() && $slice_count < $max_top_slice ) {
    $variable_start = $variable_end+1;
    $variable_end = $slice->length();

    $repeat_size = $rr->overlap_size( "1", $variable_start, $variable_end );
    my $percentage_repeat = $repeat_size / ($variable_end - $variable_start+1) * 100;

    $dfa->store( Bio::EnsEMBL::DensityFeature->new
		 (-seq_region    => $slice,
		  -start         => $variable_start,
		  -end           => $variable_end,
		  -density_type  => $variable_density_type,
		  -density_value => $percentage_repeat));
  }

  print STDOUT " DONE.\n";
}


sub register {
  my ($rr, $slice, $start, $end ) = @_;

  my $subSlice = $slice->sub_Slice( $start, $end );
  my $repeats = $subSlice->get_all_RepeatFeatures();
  for my $repeat ( @$repeats ) {
    $rr->check_and_register( "1", $repeat->seq_region_start(), $repeat->seq_region_end() );
  }
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
#  my $avg = undef;
#  $avg = $sum/$length if($length < 0);
#  print("Sum=$sum, Length=$length, Avg/Base=$sum");
}


sub usage {
  my $indent = ' ' x length($0);
  print <<EOF; exit(0);

What does it do?

Calculates the percentage of repetetive elements for top level seq_regions.

First it needs repeat_feature table to be populated. It then
deletes all PercentageRepeat entries from the analysis, density_type
and density_feature tables. All toplevel slices are fetched and sorted
from longest to shortest. For each toplevel slice, both the small and
the variable density repeats are calculated.

Small repeats are done like this: Move along the toplevel slice 1 KB
at a time. Find the %repeat for each of these 1 KB blocks.

Variable repeats are done like this: Divide each slice into 150
sub_slices.  Move along the toplevel slice 1 MB at a time. Foreach
sub_slice within the 1 MB, calculate the %repeat for that sub_slice.
Variable repeats are only found for the 100 longest toplevel slices.

Input data: repeat features, top level seq regions 
Output tables: updates analysis creation date,
               density_type (two entries, one for small_density type of 
               length 1 KB and one for variable_density_type of length 1MB), 
	       density_feature


When to run it in the release cycle?

It can be run after genebuilders have finished their Xrefs 
(script not affected by projected Xrefs).


Which databases to run it on?

Run on core databases for new species or if one of the following changed:
  - dna sequence
  - assembly
  - repeats


How long does it take?

It takes about 1 hour to run for a database in the long queue. The script is 
slowed down considerably for very fragmented genomes with thousands of toplevel
seqs. 


Usage:
  $0 -h host [-port port] -u user -p password \\
  $indent -d database  \\
  $indent [-mhost ensembl_production host] [-mport ensembl_production host] [-muser ensembl_production host] \\
  $indent [-help]  \\


  -h|host            Database host to connect to

  -port              Database port to connect to (default 3306)

  -u|user            Database username

  -p|pass            Password for user

  -mhost              ensembl_production database host to connect to

  -mport              ensembl_production database port to connect to

  -muser              ensembl_production database username

  -d|dbname          Database name

  -help              This message


EOF

}


  



=pod

SYNOPSIS

  Script to add StopGain or StopLoss transcript attributes for transcripts that contain a variation which causes gain or loss of a stop codon.

DESCRIPTION

  This script will add transcript attributes to transcripts in a Core database, indicating that the transcript
  contains a variation allele that causes either the loss or gain of a stop codon. The allele will have to have a
  minimum frequency and a mimimum number of observations (chromosomes) in a given population.
  
  Data can either be supplied through a tab-separated input file, or it will be fetched from a Variation database
  
EXAMPLE

  Command line example:
    perl nonsense_transcript_attribs.pl -help
    perl nonsense_transcript_attribs.pl -chost ens-genomics1 -cport 3306 -cuser ******** -cpass ********** -cdbname homo_sapiens_core_58_37c -vdbname homo_sapiens_variation_58_37c
    perl nonsense_transcript_attribs.pl -chost ens-genomics1 -cport 3306 -cuser ******** -cpass ********** -cdbname homo_sapiens_core_58_37c -vdbname homo_sapiens_variation_58_37c -min_count 20 -pop_group 1000Genomes
    
=cut

#!/usr/local/ensembl/bin/perl-w

use strict;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp); 
use DBI qw(:sql_types);
use Getopt::Long;

# A hard-coded hash for consequence type => attrib_type.code mapping
my %ATTRIBUTE_CODE = (
  'stop_lost' => 'StopLost',
  'stop_gained' => 'StopGained'
);

# The minimum frequency of a STOP_*-causing allele in a population to be included
my $FREQUENCY_CUTOFF = 0.1;
# The mimimum number of observations that the allele frequency is based on
my $OBSERVATION_COUNT_CUTOFF = 20;

# The names of the populations we will be looking at
my %POPULATIONS = (
  'HapMap'  => [
    'CSHL-HAPMAP:HAPMAP-ASW',
    'CSHL-HAPMAP:HapMap-CEU',
    'CSHL-HAPMAP:HAPMAP-CHB',
    'CSHL-HAPMAP:HAPMAP-CHD',
    'CSHL-HAPMAP:HAPMAP-GIH',
    'CSHL-HAPMAP:HapMap-HCB',
    'CSHL-HAPMAP:HapMap-JPT',
    'CSHL-HAPMAP:HAPMAP-LWK',
    'CSHL-HAPMAP:HAPMAP-MEX',
    'CSHL-HAPMAP:HAPMAP-MKK',
    'CSHL-HAPMAP:HAPMAP-TSI',
    'CSHL-HAPMAP:HapMap-YRI'
  ],
  '1000Genomes' => [
    '1000GENOMES:pilot_3_CEU_exon_capture_panel',
    '1000GENOMES:pilot_3_CHB_exon_capture_panel',
    '1000GENOMES:pilot_3_CHD_exon_capture_panel',
    '1000GENOMES:pilot_3_JPT_exon_capture_panel',
    '1000GENOMES:pilot_3_LWK_exon_capture_panel',
    '1000GENOMES:pilot_3_TSI_exon_capture_panel',
    '1000GENOMES:pilot_3_YRI_exon_capture_panel',
    '1000GENOMES:pilot_1_CEU_low_coverage_panel',
    '1000GENOMES:pilot_1_CHB+JPT_low_coverage_panel',
    '1000GENOMES:pilot_1_YRI_low_coverage_panel'
  ]
); 

# initialize variables
my $chost;
my $cuser;
my $cpass;
my $cport;
my $cdbname;
my $vhost;
my $vuser;
my $vpass;
my $vport;
my $vdbname;
my @population_groups;

my $path = 'GRCh37';
my $file;

$| = 1; # buffer
my $store;
my $clean;
my $list;
my $help;

usage() if (!scalar(@ARGV));

# options to be read in from the commandline
&GetOptions(
  'chost:s'     => \$chost,
  'cuser:s'     => \$cuser,
  'cdbname:s'   => \$cdbname,
  'cpass:s'     => \$cpass,
  'cport:n'     => \$cport,
  'vhost:s'     => \$vhost,
  'vuser:s'     => \$vuser,
  'vdbname:s'   => \$vdbname,
  'vpass:s'     => \$vpass,
  'vport:n'     => \$vport,
  'path:s'      => \$path,
  'file:s'      => \$file,
  'store'       => \$store,
  'clean!'      => \$clean,
  'min_freq=f'  => \$FREQUENCY_CUTOFF,
  'min_count=i'   => \$OBSERVATION_COUNT_CUTOFF,
  'pop_group=s@'  => \@population_groups,
  'list!'       => \$list,
  'help!'       => \$help
);

usage() if (defined($help));

# If a listing of population groups are desired, do that
if ($list) {
  map { printf("\%s\n",$_); map { printf("\t\%s\n",$_) } @{$POPULATIONS{$_}}; } keys(%POPULATIONS);
  exit;
}

# some checks, if not explicitly specified otherwise, will use the same connection credentials for variation as for core but the database name must be specified
if (!defined $chost || !defined $cuser || !defined $cdbname || !defined $cpass || !defined $cport || (!defined($file) && !defined $vdbname)) {
  throw("Please enter -chost -cuser -cdbname -cpass -cport -vdbname");
}
# If population group has been specified but does not match the available ones, throw exception
if (@population_groups && grep {!exists($POPULATIONS{$_})} @population_groups) {
  throw("The entered population group was not recognized. Possible choices are: " . join(",",keys(%POPULATIONS)));
}

# If no population groups were specified, use all available
unless (@population_groups) {
  @population_groups = keys(%POPULATIONS);
}

# Use the same db credentials for variation as for core if not specified otherwise
$vhost ||= $chost;
$vuser ||= $cuser;
$vpass = $cpass if (!defined($vpass));
$vport ||= $cport;

# connect to core db
# assume it has dna in it
my $cdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                            -host => $chost,
                                            -user => $cuser,
                                            -pass => $cpass,
                                            -port => $cport,
                                            -dbname => $cdbname
) or die("Could not get a database adaptor to $cdbname on $chost");
print "Connected to " . $cdb->dbc->dbname() . " on " . $cdb->dbc->host() . ":" . $cdb->dbc->port() . "\n";

# connect to the variation db if no input file has been specified
my $vdb;
if (!defined($file)) {
  $vdb = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
                                            -host => $vhost,
                                            -user => $vuser,
                                            -pass => $vpass,
                                            -port => $vport,
                                            -dbname => $vdbname
  ) or die("Could not get a database adaptor to $vdbname on $vhost");
  print "Connected to " . $vdb->dbc->dbname() . " on " . $vdb->dbc->host() . ":" . $vdb->dbc->port() . "\n";
  
  # Check that the Core and Variation databases are on the same release, throw an error if they are not
  my $cversion = $cdb->get_MetaContainerAdaptor()->get_schema_version();
  my $vversion = $vdb->get_MetaContainerAdaptor()->get_schema_version();
  throw("Core database is release $cversion but Variation database is release $vversion. The releases should match") if ($cversion != $vversion);
}

# If -clean has not been specified but -store has, check that the database does not already contain old transcript_attribs for STOP_*. These should be cleaned out.
if (!$clean && $store) {
  
  # Get a condition string for the attribute codes
  my $attrib_codes = "at.code = '" . join("' OR at.code = '",values(%ATTRIBUTE_CODE)) . "'";
  
  # Count the number of transcripts with STOP_* attributes assigned. This should always be 0 
  my $stmt = qq{
    SELECT
      COUNT(DISTINCT ta.transcript_id)
    FROM
      attrib_type at JOIN transcript_attrib ta USING (attrib_type_id)
    WHERE
      $attrib_codes
  };
  my $transcripts = $cdb->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
  
  # If the query failed, do not proceed
  die("Could not determine if transcript attributes already exist in database. You need to debug the script!") if (!defined($transcripts));
  
  # If transcripts with attributes exist, throw an error 
  throw("You have not specified old transcript attributes to be cleared but such attributes exist for $transcripts transcripts in the database. Please re-run the script with the -clean option in order to first delete those attributes") if ($transcripts > 0);
}

# If -clean has been specified, remove all transcript_attrib entries with attrib_type_id corresponding to STOP_* events
if ($clean) {
  
  # Get a condition string for the attribute codes
  my $attrib_codes = "at.code = '" . join("' OR at.code = '",values(%ATTRIBUTE_CODE)) . "'";
  
  my $stmt = qq{
    DELETE FROM
      ta
    USING
      attrib_type at JOIN transcript_attrib ta USING (attrib_type_id)
    WHERE
      $attrib_codes
  };
  my $result = $cdb->dbc->do($stmt);
  die("Failed to delete transcript_attrib entries!") if (!defined($result));
}

# Put the null transcripts into an array where each tab-separated element has the form [Transcript id] [Population] [rs-id] [Consequence]
# ENST00000492164 CSHL-HAPMAP:HapMap-YRI  rs1131265       STOP_LOST
# ENST00000493812 CSHL-HAPMAP:HapMap-HCB  rs2273865       STOP_GAINED
# ENST00000493812 CSHL-HAPMAP:HapMap-JPT  rs2273865       STOP_GAINED
# ENST00000493812 CSHL-HAPMAP:HapMap-YRI  rs2273865       STOP_GAINED
my @null_transcripts;

# If no input file has been specified, get the data from the variation database
if (!defined($file)) {
  
    # Get the source_id for dbSNP
    my $stmt = qq{
        SELECT
            source_id
        FROM
            source
        WHERE
            name = 'dbSNP'
        LIMIT 1
    };
    my $source_id = $vdb->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    # Get the sample_ids for the populations
    my %samples;
    $stmt = qq{
        SELECT
            p.sample_id,
            s.name
        FROM
            sample s JOIN
            population p ON (
                p.sample_id = s.sample_id
            )
        WHERE
            s.name IN 
    };
    foreach my $pop_group (@population_groups) {
      my $pop_stmt = $stmt . " ('" . join("','",@{$POPULATIONS{$pop_group}}) . "')";
      map {$samples{$_->[0]} = $_->[1]} @{$vdb->dbc->db_handle->selectall_arrayref($pop_stmt)};
    }
    my $sample_ids = join(",",keys(%samples));
  
    # A prepared statement for getting the populations where the stop_* causing allele has enough frequency
    $stmt = qq{
      SELECT
        a.sample_id
      FROM
        allele a JOIN
        variation v ON (
            v.variation_id = a.variation_id
        )
      WHERE
        a.variation_id = ? AND
        a.allele = ? AND
        a.frequency >= $FREQUENCY_CUTOFF AND
        a.count >= $OBSERVATION_COUNT_CUTOFF AND
        a.sample_id IN ($sample_ids) AND
        v.source_id = $source_id
    };
    my $pop_sth = $vdb->dbc->prepare($stmt);
    
  # Query the MySQL database directly for transcripts with consequence types STOP_*
  $stmt = qq{
    SELECT
      tv.feature_stable_id,
      vf.variation_id,
      vf.variation_name,
      tv.allele_string,
      tv.consequence_types
    FROM
      transcript_variation tv JOIN
      variation_feature vf USING (variation_feature_id)
    WHERE
        (
            FIND_IN_SET('stop_lost',tv.consequence_types) OR
            FIND_IN_SET('stop_gained',tv.consequence_types)
        ) AND
        tv.somatic = 0
  };
  my $sth = $vdb->dbc->prepare($stmt) or die("Error preparing statement $stmt");
  $sth->execute();
  
  # For each variation_feature, check that the source, population and minor allele frequencies are 'dbSNP', 'CSHL-HapMap-*' and '>= $FREQUENCY_CUTOFF'
  while (my ($transcript_stable_id,$variation_id,$variation_name,$allele_string,$consequence_string) = $sth->fetchrow_array()) {
    
    # Get the allele and the consequence
    my ($allele) = $allele_string =~ m/\/(.*)$/;
    my ($consequence) = $consequence_string =~ m/(stop_[^,]+)/i;
    
    # Get the HapMap population and stop_* event causing allele that has an allele frequency above the $FREQUENCY_CUTOFF
    $pop_sth->execute($variation_id,$allele);
    
    # For each result, push a tab-separated result string into the null_transcripts-array
    while (my ($sample_id) = $pop_sth->fetchrow_array()) {
      my $str = join("\t",($transcript_stable_id,$samples{$sample_id},$variation_name,$consequence));
      push(@null_transcripts,$str);
    }
    
  }
  
  $sth->finish();
}
# Otherwise, read entries from the input file and put them in the array
else {
  open(INFILE, "<$file") or die ("Can't read $file $! \n");
  LINE: while (<INFILE>) {
    my $line = $_;
    chomp $line;
  
    # ignore comments
    if ($line =~ /^#/){
      next;
    }
    
    push(@null_transcripts,$line);
  }
}

my %unique;

# Process each element in the array
while (my $line = shift(@null_transcripts)) {
  my $transc_stable_id;
  my $rsid;
  my $population;
  my $consequence;

  ($transc_stable_id,$population,$rsid,$consequence) = $line =~ m/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$/;
  print join("\t",($transc_stable_id,$population,$rsid,$consequence)) . "\n";
  
  if (!$transc_stable_id || !$rsid || !$population || !$consequence) {
    throw("Have not found all variables: $line");
  }
  print "  FOUND: $transc_stable_id $rsid $population $consequence\n";
  
  if (!$unique{$transc_stable_id.":".$rsid.":".$population.":".$consequence}) {
    store_attribute($cdb, $transc_stable_id, $rsid, $population, $consequence, $store);
    $unique{$transc_stable_id.":".$rsid.":".$population.":".$consequence} = 1;
  } else {
    print STDERR "  Already seen\n";
  }
}

sub store_attribute {
  my ($db, $stable_id, $rs_id, $pop, $cons, $do_store) = @_;

  # fetch transript
  my $transcript = $db->get_TranscriptAdaptor->fetch_by_stable_id($stable_id); 
  if (!$transcript) {
    throw("Transcript $stable_id not found");
  }
  
  # decide on which code to use
  my ($code,$name,$description) = fetch_attrib_type_by_code($db,$cons);
  if (!$code || !$name || !$description) {
    throw("Attribute type not found for cons $cons");
  }

  # now make the transcript attrib
  my $attribute = Bio::EnsEMBL::Attribute->new(
                  -CODE => $code,
                  -NAME => $name,
                  -DESCRIPTION => $description,
                  -VALUE => $rs_id.",".$pop);

  print STDERR "  Got attribute ".$attribute->code." with value ".$attribute->value." for transcript ".$transcript->stable_id."\n";

  # now store the attrib
  if ($do_store) {
    #$transcript->add_Attributes($attribute);
    #$db->get_TranscriptAdaptor->store($transcript);
    $db->get_AttributeAdaptor->store_on_Transcript( $transcript, [$attribute] );
    print STDERR "  Stored attribute ".$attribute->code." with value ".$attribute->value." on transcript ".$transcript->stable_id."\n";
  }
  return;
}

sub fetch_attrib_type_by_code {
  my ($db,$cons_code) = @_;

  my $code_tmp = $ATTRIBUTE_CODE{$cons_code};

  my $sql = "SELECT code, name, description FROM attrib_type where code = '$code_tmp'";

  my $sth = $db->dbc->prepare($sql) or die "sql error";
  $sth->execute();
  my ($code,$name,$description) = $sth->fetchrow_array;
  $sth->finish;

  return ($code,$name,$description);
}

sub usage {
	
  print qq{
  Usage: perl $0 [OPTION]
  
  Add transcript attributes to transcripts in a Core database, indicating that the transcript
  contains a variation allele that causes either the loss or gain of a stop codon. The allele will
  have to have a minimum frequency and observation count in a particular "population group" 
  (1000Genomes or HapMap).
  	
  Options:
  
      -clean      If specified, transcript attributes for stop codon losses
                  or gains already present in the Core database will be
                  deleted. This is required before storing new attributes (Optional)
                
      -store      If specified, new transcript attributes will be stored in
                  the core database. This requires all old attributes to have
                  been deleted (use -clean option) (Optional)
                
      -min_freq   The minimum frequency cutoff for an allele causing a stop loss or
                  gain event to be included. Default = 0.1 (Optional)
                
      -min_count  The minimum number of observations on which the frequency calculation 
                  is based. Default = 20 (Optional)
                
      -pop_group  The analysis can be limited to particular groups of populations, e.g.
                  'HapMap' or '1000Genomes'. By default all available population groups
                  will be analyzed (Optional).
                
    Database connection parameters are specified on the command line
    
      -chost    Core database host      (Required)
      -cport    Core database port      (Required)
      -cdbname  Core database name      (Required)
      -cuser    Core database user      (Required)
      -cpass    Core database password  (Required)
    
    A tab-separated input file containing null transcript data can be
    specified. If this option is used, transcript attributes will be added
    using this data instead of getting it from a Variation database
    
      -file     Tab-separated input file (Optional)            
    
    If no input file is specified, null transcripts will be fetched
    from a Variation database. If the optional Variation database
    credentials are not specified, will use the same as for the Core
    database
    
      -vdbname  Variation database name     (Required)
      -vhost    Variation database host     (Optional)
      -vport    Variation database port     (Optional)
      -vuser    Variation database user     (Optional)
      -vpass    Variation database password (Optional)
      
      -list     Lists the available population groups and what populations are included
      
      -help     Print this message
  };
  
  exit(0);
}

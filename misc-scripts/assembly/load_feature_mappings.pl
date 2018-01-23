#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


=head1 NAME

load_feature_mappings.pl - load feature (genes and SNPs) location in the old
database (database with alternative assembly),  mapping from old to new assembly 
and feature location in the new database (database with reference assembly) 
to a test database
Run script test_assembly_mapping.pl to produce statistics on the quality of the 
mappings.
Sample config file load_feature_mappings.ini.example.

=head1 SYNOPSIS

load_feature_mappings.pl [arguments]

Required arguments:

  --host=HOST                 new core db host HOST
  --port=PORT                 new core db port PORT
  --user=USER                 new core db username USER
  --pass=PASS                 new core db passwort PASS
  --dbname=NAME               new core db name NAME
  --altdbname=NAME            old core db name NAME  
  --chromosomes, --chr=LIST   'all' or comma separated list of chromosomes
  --assembly=NAME             new assembly NAME
  --altassembly=NAME          old assembly NAME
  --features, --ft=LIST       features to map:genes,SNPs
  
  --vardbname                 new varation db NAME (if mapping SNPs)
  --varaltdbname              old variation db NAME (if mapping SNPs)

  if not dry run:
  --testdbname=NAME           test database name NAME (test database will be dropped and recreated)

Optional arguments:

  --conffile=filename     read parameters from FILE
                                        (default: conf/Conversion.ini)

                          (if different from --host, --port, --user, --pass):
  --althost=hOST          old core db host HOST
  --altport=PORT          old core db port PORT
  --altuser=USER          old core db username USER
  --altpass=PASS          old core db passwort PASS

  --testhost=HOST         test db host HOST
  --testport=PORT         test db port PORT
  --testuser=USER         test db username USER
  --testpass=PASS         test db passwort PASS

  --varhost=hOST          new variation db host HOST
  --varport=PORT          new variation db port PORT
  --varuser=USER          new variation db username USER
  --varpass=PASS          new variation db passwort PASS

  --varalthost=hOST       old variation db host HOST
  --varaltport=PORT       old variation db port PORT
  --varaltuser=USER       old variation db username USER
  --varaltpass=PASS       old variation db passwort PASS

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            dont write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script will map genes for a given chromosome list from the old to the new assembly and 
find gene locations in the new database using stable_ids.
It will load the results into the new database.

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

use vars qw($SERVERROOT);

sub map_slice;
sub get_old_new_feature;

BEGIN {
    $SERVERROOT = "$Bin/../";
}


$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'althost=s',
    'altport=n',
    'altuser=s',
    'altpass=s',
    'altdbname=s',
    'assembly=s',
    'altassembly=s',
    'testhost=s',
    'testport=n',
    'testuser=s',
    'testpass=s',
    'testdbname=s',
    'chromosomes|chr=s@',
    'features|ft=s@',
    'varhost=s',
    'varport=n',
    'varuser=s',
    'varpass=s',
    'vardbname=s',
    'varalthost=s',
    'varaltport=n',
    'varaltuser=s',
    'varaltpass=s',
    'varaltdbname=s',
);
$support->allowed_params(
    $support->get_common_params,
    'althost',
    'altport',
    'altuser',
    'altdbname',
    'assembly',
    'altassembly',
    'testhost',
    'testport',
    'testuser',
    'testpass',
    'testdbname',
    'chromosomes',
    'features',
    'varhost',
    'varport',
    'varuser',
    'varpass',
    'vardbname',
    'varalthost',
    'varaltport',
    'varaltuser',
    'varaltpass',
    'varaltdbname',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
  'altdbname',
  'assembly',
  'altassembly',
  'chromosomes',
  'features',
  'dbname',
  'host',
  'user'
);

if ($support->param('dry_run')) {
  $support->log_stamped("Dry run. Test database will not be updated. Writing results to logfile instead.\n\n"); 
  $support->param('verbose',1);
} else {
    if ( !$support->param('testdbname') ) {
	$support->log_error("testdbname is required if not a dry run\n",1);
    }
}

my %valid_features = ('genes' => 1,
'SNPs' => 1 );

$support->comma_to_list('features');
my @feature_types = $support->param('features');
my %feature_types = map { $_ => 1} @feature_types;

foreach my $feature_type (@feature_types ) {
    if (!exists($valid_features{$feature_type}) ) {
	delete $feature_types{$feature_type};
	print STDERR "invalid feature type: $feature_type\n";
    }
} 

if (scalar(keys %feature_types) == 0) {
    $support->log_error("listed feature(s) not valid: use genes or SNPs only\n",1);
}

if (exists($feature_types{'SNPs'}) ){
    if( !$support->param('vardbname') or !$support->param('varaltdbname') ) {
	$support->log_error("vardbname and varaltdbname are required if mapping SNPs\n",1);
    }

}


#####
# connect to database and get adaptors
#
my ($dba, $dbh, $sth);

if ( !defined($support->param('pass')) ) {
    $support->param('pass', '');
}


# first set connection parameters for alternative db and test db
if ( !defined($support->param('althost')) ) { $support->param('althost',$support->param('host')); }
if ( !defined($support->param('altport')) ) { $support->param('altport',$support->param('port')); }
if ( !defined($support->param('altuser')) ) { $support->param('altuser',$support->param('user')); }
if ( !defined($support->param('altpass')) ) { $support->param('altpass',$support->param('pass')); }

if ( !defined($support->param('testhost')) ) { $support->param('testhost',$support->param('host')); }
if ( !defined($support->param('testport')) ) { $support->param('testport',$support->param('port')); }
if ( !defined($support->param('testuser')) ) { $support->param('testuser',$support->param('user')); }
if ( !defined($support->param('testpass')) ) { $support->param('testpass',$support->param('pass')); }

# reference database
$dba->{'ref'} = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $support->param('host'),
                                            -user   => $support->param('user'),
                                            -pass   => $support->param('pass'),
                                            -port   => $support->param('port'),
                                            -dbname => $support->param('dbname')  );

$dbh->{'ref'} = $dba->{'ref'}->dbc->db_handle;

# database containing the alternative assembly
$dba->{'alt'} = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $support->param('althost'),
                                            -user   => $support->param('altuser'),
                                            -pass   => $support->param('altpass'),
                                            -port   => $support->param('altport'),
                                            -dbname => $support->param('altdbname')  );

$dbh->{'alt'} = $dba->{'alt'}->dbc->db_handle;

#test database
my $test_dbname = $support->param('testdbname');
$support->param('testdbname','');
$dbh->{'test'} = $support->get_dbconnection('test');

$support->comma_to_list('chromosomes');
my @chrs = $support->param('chromosomes');

#get chrmosomes

my %db_chrs;

my $chr_name;
$sth =  $dbh->{'ref'}->prepare("select s.name from seq_region s join coord_system c using(coord_system_id) where rank=1 and length(s.name < 3) order by s.name+0");
$sth->execute;
$sth->bind_columns(\$chr_name);
while ($sth->fetch){
  $db_chrs{$chr_name} = 1;
}
$sth->finish;

if ($chrs[0] eq 'all') {
  @chrs = sort keys %db_chrs;
} else {
  my %temp_chrs;
  foreach my $chr (@chrs) {
      if ( exists($db_chrs{$chr}) ) {
          $temp_chrs{$chr} = 1;
      } else {
        print STDERR "unknown chromosome $chr\n";
      }
  }
  @chrs = sort keys %temp_chrs;
}

if ( scalar(@chrs) == 0) {
  $support->log_error("listed chromosome(s) not found in the database\n",1);
}


if (!$support->param('dry_run')) {
    #create the test database
   
    eval{
	$dbh->{'test'}->do("drop database if exists $test_dbname");
	$dbh->{'test'}->do("create database $test_dbname");
	$dbh->{'test'}->do("use $test_dbname");
	$dbh->{'test'}->do("CREATE TABLE mapping (
  feature_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  feature_type char(30) DEFAULT NULL,
  stable_id varchar(128) DEFAULT NULL,
  old_assembly varchar(30) DEFAULT NULL,
  old_chr char(3) DEFAULT NULL,
  old_start int(20) DEFAULT NULL,
  old_end int(20) DEFAULT NULL,
  old_length int(20) DEFAULT NULL,
  old_strand tinyint(2) DEFAULT NULL,
  feature_found_in_new_db tinyint(1) unsigned DEFAULT NULL,
  new_chr char(3) DEFAULT NULL,
  new_start int(20) DEFAULT NULL,
  new_end int(20) DEFAULT NULL,
  new_length int(20) DEFAULT NULL,
  new_strand tinyint(2) DEFAULT NULL,
  new_assembly varchar(30) DEFAULT NULL,
  mapping_quality tinyint(1) unsigned DEFAULT NULL,
  mapping_start int(20) DEFAULT NULL,
  mapping_end int(20) DEFAULT NULL,
  mapping_length int(20) DEFAULT NULL,
  mapping_strands varchar(5) DEFAULT NULL,
  mapping_chrs varchar(10) DEFAULT NULL,
  PRIMARY KEY (feature_id))");
 
    };
    if ( $@ ) {
	$support->log_error("errors encountered when creating the test database\n",1);
    }
}


if ( exists($feature_types{'SNPs'}) ) {
    #connect to the variation db

    if ( !defined($support->param('varhost')) ) { $support->param('varhost',$support->param('host')); }
    if ( !defined($support->param('varport')) ) { $support->param('varport',$support->param('port')); }
    if ( !defined($support->param('varuser')) ) { $support->param('varuser',$support->param('user')); }
    if ( !defined($support->param('varpass')) ) { $support->param('varpass',$support->param('pass')); }

    if ( !defined($support->param('varalthost')) ) { $support->param('varalthost',$support->param('host')); }
    if ( !defined($support->param('varaltport')) ) { $support->param('varaltport',$support->param('port')); }
    if ( !defined($support->param('varaltuser')) ) { $support->param('varaltuser',$support->param('user')); }
    if ( !defined($support->param('varaltpass')) ) { $support->param('varaltpass',$support->param('pass')); }


						  
    $dba->{'var'} = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
	-host => $support->param('varhost'),
	-port => $support->param('varport'),
	-user => $support->param('varuser'),
	-pass => $support->param('varpass'),
	-dbname => $support->param('vardbname'));

    $dba->{'varalt'} = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
	-host => $support->param('varalthost'),
	-port => $support->param('varaltport'),
	-user => $support->param('varaltuser'),
	-pass => $support->param('varaltpass'),
	-dbname => $support->param('varaltdbname'));

}


#connect to the test mapping db to store results
if ( !$support->param('dry_run') ) {
    $sth = $dbh->{'test'}->prepare("INSERT INTO mapping(feature_type, stable_id, old_assembly,old_chr,old_start,old_end,old_length,old_strand,feature_found_in_new_db,new_chr,new_start,new_end,new_length,new_strand,new_assembly,mapping_quality,mapping_start,mapping_end,mapping_length,mapping_strands,mapping_chrs ) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");

}

my $new_sa = $dba->{'ref'}->get_adaptor("slice");
my $old_sa = $dba->{'alt'}->get_adaptor("slice");

if ( exists($feature_types{'genes'}) ) {

 
    my $new_ga = $dba->{'ref'}->get_adaptor("gene");    
    my $old_ga = $dba->{'alt'}->get_adaptor("gene");

    foreach my $chr (@chrs) {
	my $feature_type = 'gene';
	get_old_new_feature($old_ga,$new_ga,$feature_type,$chr, undef);
 
    }

}

if ( exists($feature_types{'SNPs'}) ) {

    my $new_vfa = $dba->{'var'}->get_adaptor('variationfeature');
    my $new_va =  $dba->{'var'}->get_adaptor('variation');
    my $old_vfa = $dba->{'varalt'}->get_adaptor('variationfeature');

    foreach my $chr (@chrs) {
	my $feature_type = 'SNP';
	get_old_new_feature($old_vfa,$new_vfa,$feature_type,$chr, $new_va); 
    }
}

if ( !$support->param('dry_run') ){ 
    $sth->finish();
}

# finish logfile
$support->finish_log;


sub get_old_new_feature {

    my $old_adaptor = shift;
    my $new_adaptor = shift;
    my $feature_type = shift;
    my $chr = shift;
    my $new_va = shift;

    $support->log_stamped("fetching $feature_type mappings for chromosome $chr\n",1);

    #get all genes on the chromosome from the old database
    my $old_slice = $old_sa->fetch_by_region('chromosome', $chr, undef, undef,
    undef, $support->param('altassembly'));
  
    # get a slice on the old assembly from the new database
    my $newdb_oldasm_chr = $new_sa->fetch_by_region('chromosome', $chr, undef, undef,
    undef, $support->param('altassembly'));
  
    my @old_features = @{$old_adaptor->fetch_all_by_Slice($old_slice)};
    foreach my $old_feature (@old_features) {
	my $feature_name;
	if ($feature_type eq 'gene') {
	    $feature_name = $old_feature->stable_id;
	} elsif ($feature_type eq 'SNP') {
	    $feature_name = $old_feature->variation_name;
        } 
	$support->log_verbose("$feature_type $feature_name location in old db: " . $old_feature->start . ":" . $old_feature->end . ":" . $old_feature->strand . " chromosome $chr\n",1);

	
	if ($old_feature->length == 0) {
	    next;
	}

      my $newdb_old_asm_slice = Bio::EnsEMBL::Slice->new(-coord_system => $newdb_oldasm_chr->coord_system,
                                                 -start => $old_feature->start,
                                                 -end => $old_feature->end,
                                                 -strand => $old_feature->strand,
                                                 -seq_region_name => $chr,
                                                 -seq_region_length => $old_feature->length,
                                                 -adaptor => $new_sa);
      
      my $mapping_start = $newdb_oldasm_chr->end;
      my ($mapping_end, $mapping_length, $mapping_strands, $mapping_chrs, $mapping);
      ($mapping_start, $mapping_end, $mapping_length, $mapping_strands,$mapping_chrs, $mapping) = map_slice($newdb_old_asm_slice, $mapping_start);

      #get the feature from the new db
	my $new_feature;
	if ($feature_type eq 'gene') {
	    $new_feature = $new_adaptor->fetch_by_stable_id($feature_name);
	} elsif ($feature_type eq 'SNP') {
	    
	    my $variation = $new_va->fetch_by_name($feature_name); 
	    if ($variation) {
		my @new_features = @{$new_adaptor->fetch_all_by_Variation($variation)};
		foreach my $vf (@new_features) {
		    if ($vf->allele_string eq $old_feature->allele_string) {
			$new_feature = $vf;
		    }
		}
	    }
        } 


      my $new_feature_found = 0;
      my ($new_chr,$new_start,$new_end,$new_length,$new_strand);
      if (defined($new_feature) ) {
         $new_chr = $new_feature->slice->seq_region_name;
         $new_start = $new_feature->start;
         $new_end = $new_feature->end;
         $new_strand = $new_feature->strand;
	 $new_length =  $new_feature->end - $new_feature->start + 1;
         $new_feature_found = 1;
         $support->log_verbose("$feature_type $feature_name found in the new db; location:$new_start:$new_end:$new_strand chromosome $new_chr\n",1);
 
      } else {
         $support->log_verbose("$feature_type $feature_name not found in the new db; can't compare $feature_type location to mapped location\n",1);
      }

      #store mapping in the db  
      if ( !$support->param('dry_run') ) { 
	  $sth->execute($feature_type,$feature_name,$support->param('altassembly'),$chr,$old_feature->start, $old_feature->end, $old_feature->end - $old_feature->start + 1, $old_feature->strand, $new_feature_found, $new_chr, $new_start, $new_end,$new_length, $new_strand, $support->param('assembly'), $mapping, $mapping_start, $mapping_end, $mapping_length, $mapping_strands, $mapping_chrs );
      }
 
   }

}


sub map_slice {

      my $newdb_old_asm_slice = shift;
      my $new_asm_start = shift;

      $support->log_verbose("projection to new assembly: \n",1);
      #project to new assembly
      my @segments = @{ $newdb_old_asm_slice->project('chromosome', $support->param('assembly')) };
 
      my $new_asm_end = 0;
      my %mapping_strands;
      my %mapping_chrs;
      my $first_strand;
      my $first_chr;
      my $mapping_length = 0;

      foreach my $segment (@segments) {
        my $p_slice = $segment->to_Slice;
	if (! defined($first_strand)) {
	    $first_strand = $p_slice->strand;
	}
	if (! defined($first_chr) ){
	    $first_chr = $p_slice->seq_region_name;
	}
        if ($p_slice->start < $new_asm_start && $p_slice->strand == $first_strand && $p_slice->seq_region_name eq $first_chr) {
          $new_asm_start = $p_slice->start;
        }
        if ($p_slice->end > $new_asm_end && $p_slice->strand == $first_strand && $p_slice->seq_region_name eq $first_chr ) {
          $new_asm_end =  $p_slice->end;
        }
	if ($new_asm_end > 0  && $p_slice->strand == $first_strand && $p_slice->seq_region_name eq $first_chr ) {

	    $mapping_length += $p_slice->length;
	}
        $mapping_strands{$p_slice->strand} = 1;
        $mapping_chrs{$p_slice->seq_region_name} = 1;
        $support->log_verbose($p_slice->start() . ":" . $p_slice->end() . ":" . $p_slice->strand . " chromosome ".$p_slice->seq_region_name ."\n",1);    
      }
            
      my $mapping; #0-no mapping, 1-gaps exist, 2-no gaps     

      if (scalar(@segments) == 0) {
        $support->log_verbose("no mapping available\n",1);
        $mapping = 0;
	undef $new_asm_start;
	undef $new_asm_end;
	undef $mapping_length;
      } else {   
	 $support->log_verbose("maping start $new_asm_start mapping end $new_asm_end\n",1);
	 
        if (scalar(@segments) > 1) {
          $support->log_verbose("gaps in new assembly sequence mapping\n",1);
          $mapping = 1;
        } else {
          $support->log_verbose("no gaps in new assembly sequence mapping\n",1);
          $mapping = 2;
        }
      }

      my $mapping_strands;
      if (scalar(keys %mapping_strands) > 1) {
        $mapping_strands = 'both';
      }
      if (scalar(keys %mapping_strands) == 1) {
	my @strands = keys %mapping_strands;
        $mapping_strands = $strands[0];
      }
      my $mapping_chrs;
      if (scalar(keys %mapping_chrs) >= 1) {
        	my @chrs = keys %mapping_chrs;
        $mapping_chrs = join(',',@chrs);
      }

      return($new_asm_start, $new_asm_end, $mapping_length, $mapping_strands,$mapping_chrs, $mapping);
}

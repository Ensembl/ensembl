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


=head1

test_assembly_mapping.pl [arguments]

Produce statistics on the quality of feature mappings between two different assemblies.
Feature mappings need to be loaded first using script load_feature_mappings.pl
Sample config file: test_assembly_mapping.ini.example

Required arguments:
  --host=HOST                 test database host HOST
  --port=PORT                 test database port PORT
  --user=USER                 test database username USER
  --pass=PASS                 test database passwort PASS
  --dbname=NAME               test database name NAME


Optional arguments:

  --conffile=filename     read parameters from FILE
                                        (default: conf/Conversion.ini)
  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../";
}


$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);

$support->allowed_params(qw(dbname host port user pass conffile logfile logpath logappend));

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(qw(dbname host port user));

my ($dba, $dbh, $sql, $sth);
$dbh = $support->get_dbconnection();


my %mapping_quality = ( 1 => 'with gaps',
			2 => 'without gaps');

#get chromosomes

my @chrs = @{$dbh->selectcol_arrayref("select distinct old_chr from mapping order by old_chr+0")};

foreach my $chr (@chrs) {
 
  my @feature_types = @{$dbh->selectcol_arrayref("select distinct feature_type from mapping where old_chr = '". $chr . "' order by feature_type")};

  foreach my $feature_type (@feature_types) {
      $support->log_stamped("processing $feature_type mappings for chromosome $chr\n\n",1);


      #count all features in the old db
      $sql = "select count(1) from mapping where feature_type = '" . $feature_type . "' and old_chr = '". $chr . "'" ;
      my ($total_features_old) = $dbh->selectrow_array($sql);

      #count all features both in the old and new db
      $sql = "select count(1) from mapping where feature_type = '" . $feature_type . "' and feature_found_in_new_db = 1 and old_chr = '". $chr . "'";
      my ($total_features_old_new) = $dbh->selectrow_array($sql);

      my $perc_features = ($total_features_old_new/$total_features_old) * 100;
      $perc_features = sprintf "%.2f", $perc_features;

      $support->log_stamped("$total_features_old_new $feature_type" . "s found in both old and new assembly dbs ($perc_features\% in old assembly db found in new db)\n\n", 1);

      #count features which have the same length in the new db
      $sql = "select count(1) from mapping where feature_type = '" . $feature_type . "' and  feature_found_in_new_db = 1 and old_length = new_length and old_chr = '". $chr . "'";

      my ($total_same_length) = $dbh->selectrow_array($sql);

      my $ok_mappings_1 = 0;
      my $perc_mappings_1 = 0;

      if ($total_same_length > 0) {
	  #count mappings where start and end match feature location in the new db by mapping quality (1-contains gaps, 2 - no gaps

	  $support->log_stamped("$total_same_length $feature_type" . "s where feature length in the old assembly db is the same as in the new assembly db\n\n", 2);

	  $sql = "select count(1) as mapping_count,mapping_quality from mapping where feature_type = '" . $feature_type . "' and old_chr = '". $chr . "' and  feature_found_in_new_db = 1 and old_length = new_length and new_start = mapping_start and new_end = mapping_end and (new_strand = mapping_strands or mapping_strands = 'both') and (mapping_chrs = new_chr or mapping_chrs like concat(new_chr,',%')) group by mapping_quality";

	  my %mapping_counts = %{$dbh->selectall_hashref($sql,"mapping_quality")};
	  
	  my $ok_mappings = 0;

	  foreach my $mapping_q (keys %mapping_counts) {
		   my $perc = ($mapping_counts{$mapping_q}{'mapping_count'}/$total_same_length) * 100;
		   $perc = sprintf "%.2f", $perc;

		   $support->log_stamped("$mapping_counts{$mapping_q}{'mapping_count'} $feature_type mappings $mapping_quality{$mapping_q} where feature length in the old assembly db is the same as in the new assembly db and mapping start and end match feature location in the new db ($perc\%)\n\n", 3);
		  
		   $ok_mappings += $mapping_counts{$mapping_q}{'mapping_count'};
		   
	  }

	  $ok_mappings_1 += $ok_mappings;
	  my $perc_mappings = ($ok_mappings/$total_same_length) * 100;
	  $perc_mappings = sprintf "%.2f", $perc_mappings;

	  $support->log_stamped("total: $ok_mappings ok $feature_type mappings ($perc_mappings\%):\n\n", 3);

      } else {
	  $support->log_stamped("no $feature_type" . "s found for chromosome $chr where feature length in the old assembly db is the same as in the new assembly db\n\n", 2);
      }
   
      $sql = "select count(1) from mapping where feature_type = '" . $feature_type . "' and old_chr = '". $chr . "' and feature_found_in_new_db = 1 and old_length != new_length";

      my ($total_diff_length) = $dbh->selectrow_array($sql);

      if ($total_diff_length > 0) {      

	  $support->log_stamped("$total_diff_length $feature_type" . "s found for chromosome $chr where feature length in the old assembly db is different than in the new assembly db\n\n", 2);

	  #count features with different length in the new db where either mapping start or end match new feature start or end respectively, by mapping quality
	  
	  $sql = "select count(1) as mapping_count,mapping_quality from mapping where feature_type = '" . $feature_type . "' and old_chr = '". $chr . "'  and  feature_found_in_new_db = 1 and old_length != new_length and (new_start = mapping_start or new_end = mapping_end) and (new_strand = mapping_strands or mapping_strands = 'both') and (mapping_chrs = new_chr or mapping_chrs like concat(new_chr,',%')) group by mapping_quality";

	  my %mapping_counts = %{$dbh->selectall_hashref($sql,"mapping_quality")};

	  my $ok_mappings = 0;
	 
	  foreach my $mapping_q (keys %mapping_counts) {
	       my $perc = ($mapping_counts{$mapping_q}{'mapping_count'}/$total_diff_length) * 100;
	       $perc = sprintf "%.2f", $perc;

	       $support->log_stamped("$mapping_counts{$mapping_q}{'mapping_count'} $feature_type mappings $mapping_quality{$mapping_q} found for chromosome where feature length in the old assembly db is different than in the new assembly db and either mapping start or end matches the gene location in the new db ($perc\%)\n\n", 3);

	       $ok_mappings += $mapping_counts{$mapping_q}{'mapping_count'};
	       
	  }

	  #count features with different length in the new db where either mapping start or end is within the difference between old and new length from the new feature start or end respectively, by mapping quality
	  $sql = "select count(1) as mapping_count,mapping_quality from mapping where feature_type = '" . $feature_type . "' and old_chr = '". $chr . "' and  feature_found_in_new_db = 1 and old_length != new_length and (abs(new_start - mapping_start) <= abs(new_length - old_length) or abs(new_end - mapping_end) <= abs(new_length - old_length)) and new_start != mapping_start and new_end != mapping_end and (new_strand = mapping_strands or mapping_strands = 'both') and (mapping_chrs = new_chr or mapping_chrs like concat(new_chr,',%')) group by mapping_quality";

	  %mapping_counts = %{$dbh->selectall_hashref($sql,"mapping_quality")};


	  foreach my $mapping_q (keys %mapping_counts) {
	       my $perc = ($mapping_counts{$mapping_q}{'mapping_count'}/$total_diff_length) * 100;
	       $perc = sprintf "%.2f", $perc;

	       $support->log_stamped("$mapping_counts{$mapping_q}{'mapping_count'} $feature_type mappings $mapping_quality{$mapping_q} found for chromosome where feature length in the old assembly db is different than in the new assembly db and either mapping start or end is within the difference between old and new feature length from the new feature start or end, respectively ($perc\%)\n\n", 3);
	       $ok_mappings += $mapping_counts{$mapping_q}{'mapping_count'};
	   
	  }

	  my $perc_mappings = ($ok_mappings/$total_diff_length) * 100;
	  $perc_mappings = sprintf "%.2f", $perc_mappings;
	  $support->log_stamped("total: $ok_mappings ok $feature_type mappings ($perc_mappings\%):\n\n", 3);

	  $ok_mappings_1 += $ok_mappings;	  

      }	 else {
	  $support->log_stamped("no $feature_type" . "s found for chromosome $chr where feature length in the old assembly db is different than in the new assembly db\n\n", 2);
      } 

      $perc_mappings_1 = ($ok_mappings_1/$total_features_old_new) * 100;
      $perc_mappings_1 = sprintf "%.2f", $perc_mappings_1;
      $support->log_stamped("total: $ok_mappings_1 ok $feature_type mappings ($perc_mappings_1\%):\n\n", 1);
      

  }
 
}

# finish logfile
$support->finish_log;

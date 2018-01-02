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

use strict;
use warnings;

use Getopt::Long;

my $path;
my $host;
my $user;

my $ret = Getopt::Long::GetOptions ('path=s' => \$path,
                                    'host=s' => \$host,
                                    'user=s' => \$user,
                                    'help'   => sub { usage(); exit(0); } );

if(!defined $path){
  print "you must defined the path to the ensembl webcode\n";
  usage();
  exit(0);
}

my $file = $path."/ensembl-webcode/htdocs/info/docs/api/core/core_tutorial.html";

my $code_count = 0;
$user ||= "ensro";
$host ||= "ens-staging";
my $header =(<<"HED");
####### start of insertion #########
use Bio::EnsEMBL::Registry;
my \$registry = 'Bio::EnsEMBL::Registry';

\$registry->load_registry_from_db(
    -host => "$host",
    -user => "$user");
HED

my $slice_adaptor = (<<'SLI');
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
SLI

open(HTML,"<$file") || die "Could not open $file\n";

my $code_mode = 0;
my $code = "";
while (my $line = <HTML>){
  chomp $line;
#  print $line."DUDE\n";
  if($code_mode){
    if($line =~ m/\<\/pre\>/){
      $code_mode = 0;
      process_code();
      $code = "";
    }
    else{
      $code .= "\n".$line;
    }
  }
  elsif($line =~ m/pre.*code/){
#    print "CODEMODE*******\n"; 
    $code_mode = 1;
  }

}

close HTML;

sub process_code{

  $code_count++;
  $file = "test".$code_count.".pl";
  open(PERL,">".$file) or die "Could not open $file for writing\n";

  $code =~ s/\&gt\;/\>/g;
  my $prefix = "";
  if(!($code =~ /load_registry_from_db/m)){
    $prefix = $header;
    if(!($code =~ /my\s*\$slice_adaptor/)){
      $prefix .= $slice_adaptor;
      if(!($code =~ /fetch_by_region/)){
	$prefix .= 'my $slice = $slice_adaptor->fetch_by_region( "clone", "AL031658.11" );'."\n";
      }
    }
    if(!($code =~ /my\s*\$cs_adaptor/)){
      $prefix .= 'my $cs_adaptor = $registry->get_adaptor( "Human", "Core", "CoordSystem" );'."\n";;
    }
    if(!($code =~ /my\s*\$gene_adaptor/)){
      $prefix .= 'my $gene_adaptor = $registry->get_adaptor( "Human", "Core", "Gene" );'."\n";;
    }
    if($code =~ /\$feature\-\>/ and !( $code =~ /\$feature_adaptor/) ){
      $prefix .= 'my $feat_adaptor = $registry->get_adaptor( "Human", "Core", "Gene" );'."\n";
      $prefix .= 'my $feature = $feat_adaptor->fetch_by_display_label("COG6");'."\n";
    }
    if($code =~ /\$transcript\-\>/ and ! ( $code =~ /\$transcript_adaptor/) ){
      $prefix .= 'my $transcript_adaptor = $registry->get_adaptor( "Human", "Core", "transcript" );'."\n";
      $prefix .= 'my $transcript = $transcript_adaptor->fetch_by_stable_id("ENST00000380152");'."\n";
    }
    if($code =~ /\$translation\-\>/ and ! ( $code =~ /\$translation_adaptor/) ){
      $prefix .= 'my $transcript_adaptor = $registry->get_adaptor( "Human", "Core", "transcript" );'."\n";
      $prefix .= 'my $transcript = $transcript_adaptor->fetch_by_stable_id("ENST00000380152");'."\n";
      $prefix .= 'my $translation = $transcript->translation;'."\n";
    }
    if( $code =~ /\$marker\-\>/ and ! ( $code =~ /\$marker_adaptor/) ){
      $prefix .= 'my $marker_adaptor = $registry->get_adaptor( "Human", "Core", "marker" );'."\n";
      $prefix .= 'my $marker = $marker_adaptor->fetch_all_by_synonym("D9S1038E")->[0];'."\n";
    }
  }
  $code = $prefix."####### end of insertion#########\n".$code;
  print PERL $code."\n";;
  close PERL;
  sleep 2;
#  print "CODE".$code."\n";

  my $test_val = system("perl $file >& /dev/null");
#  print "VAL is $test_val\n";


 # my $test_val =  eval {$code};

  if($test_val){
    print "file $file FAILED  $test_val\n";
  }
  else{
    print "file $file OKAY    $test_val\n";
  }
}

sub usage {
  print << "EOF";
  perl test_code_from_html.pl

  This script will strip and create seperate perl files for any block in the html
  that is in surrounded by <pre class="code"> ... </pre>
  The script adds registry code etc at the start if is not a full code example.
  These perl scripts are then executed and the status of this run reports if the
  code as successfully executed.

  -path    The path to the ensembl-webcode directory (Must be set)

  -user    User name to be used for the mysql instance

  -host    Host name of the mysql instance (defualt ens-staging)

  i.e. perl test_code_from_html -path ~/src

EOF
}

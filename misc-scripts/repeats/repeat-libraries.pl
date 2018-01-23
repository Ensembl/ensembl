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

#
# Repeat classification script
# based on js5's lite database repeat-libraries script
# 
# script  repeat-libraries.pl <UNUSED>
#
# This script is used to do run on (old) v19 databases to get the 
# repeat class from the repeat name before categorising them into types. 
# It is not used for any other purpose anymore. Repeat classification on 
# newer v32 databases is done with repeat-types.pl
#  

use strict;
use warnings;


use DBI;
use Getopt::Long;

my ( $host, $user, $pass, $port, $expression, $dbpattern, $repeatfile, $help );

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "repeatfile=s", \$repeatfile,
	    "dbpattern=s", \$dbpattern,
      "help", \$help
	  );

if($help) {
  usage();
}

if( !$host ) {
  print STDERR "-host argument is required\n";
  usage();
}

if( !$dbpattern ) {
  print STDERR "-dbpattern argument is required\n";
  usage();
}

if( !$repeatfile) {
  print STDERR "-repeatfile argument is required\n";
  usage();
}

my $dsn = "DBI:mysql:host=$host";
if( $port ) {
  $dsn .= ";port=$port";
}

my $dbh = DBI->connect( $dsn, $user, $pass );

my @dbnames = map {$_->[0] } @{ $dbh->selectall_arrayref( "show databases" ) };

my @dbs = grep {$_ =~ /$dbpattern/} @dbnames;

foreach my $db (@dbs) {
  open RFILE, $repeatfile or die("Could not open repeat file $repeatfile");

  print STDERR "$db\n";

  $dbh->do("use $db");

  print STDERR "  Clearing repeat_class\n";

  $dbh->do("update repeat_consensus set repeat_class = ''");

  print STDERR "  Reading specific repeat classes from input file\n";


  my $C=0;
  while(<RFILE>) {
    chomp;
    my($hid,$type) = split( /\t/, $_, 2);
    $dbh->do("update repeat_consensus set repeat_class = ? where repeat_name in (?,?,?)", {} , $type, $hid, substr($hid,0,15), "$hid-int" );
    $C++;
    print STDERR "$C\n" unless $C % 100;
  }

  close RFILE;

  print STDERR "  Consensifying remaining repeat classes\n";

  $dbh->do("update repeat_consensus set repeat_class = 'Simple_repeat'  where repeat_class= '' and repeat_name like '%)n'" );
  $dbh->do("update repeat_consensus set repeat_class = 'low_complexity'  where repeat_class= '' and repeat_name like '%-rich'" );
  $dbh->do("update repeat_consensus set repeat_class = 'low_complexity'  where repeat_class= '' and repeat_name like 'poly%'" );

  $dbh->do("update repeat_consensus set repeat_class = 'LTR/ERVL'  where repeat_class= '' and repeat_name like '%ERVL%' " );
  $dbh->do("update repeat_consensus set repeat_class = 'LTR/ERVL'  where repeat_class= '' and repeat_name like '%ERV16%' " );
  $dbh->do("update repeat_consensus set repeat_class = 'SINE/Alu'  where repeat_class= '' and repeat_name like 'Alu%' " );
  $dbh->do("update repeat_consensus set repeat_class = 'SINE/Alu'  where repeat_class= '' and repeat_name like '%F_AM%' " );
  $dbh->do("update repeat_consensus set repeat_class = 'LINE/L1'  where repeat_class= '' and repeat_name like 'L1%' " );
  $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER2_type'  where repeat_class= '' and repeat_name like 'Tigger%' " );
  $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER1_type'  where repeat_class= '' and repeat_name like 'Charlie%' " );
  $dbh->do("update repeat_consensus set repeat_class = 'DNA/Tc2'  where repeat_class= '' and repeat_name like 'HsTC%' " );


  $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER2_type'  where repeat_class= '' and repeat_name like 'MER46%' " );
  $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER2_type'  where repeat_class= '' and repeat_name like 'MER7%' " );
  $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER1_type'  where repeat_class= '' and repeat_name like 'MER91' " );
  $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER1_type'  where repeat_class= '' and repeat_name like 'MER58' " );
  $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER1_type'  where repeat_class= '' and repeat_name like 'MER63' " );
  $dbh->do("update repeat_consensus set repeat_class = 'Satellite/telomeric'  where repeat_class= '' and repeat_name like 'SUBTEL_%' " );

  $dbh->do("update repeat_consensus set repeat_class = 'trf'  where repeat_class = '' and repeat_name = 'trf' " );
  $dbh->do("update repeat_consensus set repeat_class = 'dust' where repeat_class = '' and repeat_name = 'dust'" );


  # $dbh->do("update repeat_consensus set repeat_class = 'LTR/ERVL'  where repeat_class= 'Unknown' and repeat_name like 'MER70%' " );
  # $dbh->do("update repeat_consensus set repeat_class = 'DNA/AcHobo'  where repeat_class= 'Unknown' and repeat_name like 'ORSL' " );

  print STDERR "  Setting repeat types\n";

  my %mappings = (
    'Low_Comp%' => 'Low complexity regions',
    'LINE%'	=> 'Type I Transposons/LINE',
    'SINE%'	=> 'Type I Transposons/SINE',
    'DNA%'	=> 'Type II Transposons',
    'LTR%'	=> 'LTRs',
    'Other%'	=> 'Other repeats',
    'Satelli%'	=> 'Satellite repeats',
    'Simple%'	=> 'Simple repeats',
    'Other%'	=> 'Other repeats',
    'Tandem%'	=> 'Tandem repeats',
    'TRF%'	=> 'Tandem repeats',
    'dust%' => 'Dust',
    'Unknown%'	=> 'Unknown',
    '%RNA'	=> 'RNA repeats',
  );
  foreach (keys %mappings) { 
    $dbh->do(qq(update repeat_consensus set repeat_type = '$mappings{$_}' where repeat_class like '$_')); 
  }

  # type all remaining repeats as unknown
  $dbh->do(qq(update repeat_consensus set repeat_type = 'Unknown' where repeat_type = ''));
  $dbh->do(qq(update repeat_consensus set repeat_type = 'Unknown' where repeat_type = null));
}

print STDERR "All done.\n";

$dbh->disconnect;


sub usage {
  print STDERR <<EOF

This program classifies the repeats stored in a core database into some
somewhat sensible categories.  It does this through a combination of a
repeat.txt file extracted from RepeatMasker repeat libraries and through
some simple pattern matching of the repeat names.

usage: perl repeat-libraries.pl  [-user <user>] [-port <port>] [-pass <pass>]
               -host <host> -dbpattern <regexp> -repeatfile <file>

example: perl repeat-libraries.pl -user ensadmin -pass secret -host ecs1g \\
                -dbpattern '^homo_sapiens_(core|vega)_20_34c$' -repeatfile repeats.txt

EOF
;
  exit;
}

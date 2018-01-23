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

# this script tries to build repeat consensus sequences
# for TRF repeats that are unique and minimal.This should reduce the amount 
# of TRF consensus objects.

# further it will replace all repeat consensi longer than 8 characters with 
# just one for the length.

use strict;
use warnings;

use DBI;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname );

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname
	  );

if( !$host || !$dbname ) {
  usage();
}

my $dsn = "DBI:mysql:host=$host;dbname=$dbname";
if( $port ) {
  $dsn .= ";port=$port";
}

# retrive all consensi of length 1..8
# build the normal consensus ( rotate, revcomp and build alphabetical minimum )

my $db = DBI->connect( $dsn, $user, $pass );

# check if the trfs have the uncompressed format

my $res = $db->selectall_arrayref( "select count(*) from repeat_consensus where repeat_consensus rlike \"\\\\(\"" );

if( $res->[0]->[0] > 0 ) {
  print STDERR "Database might alread be converted, contains \"(\" in repeat_consensus.\n";
  exit;
}


my $filename = "tmp_old_new_rcid.txt";
open ( FH, ">$filename" ) or die ( "Couldnt write $filename " );

for my $i ( 1..8 ) {

  my $remap_count = 0;
  my $all_cons_count = 0;
  my $sth = $db->prepare( "select repeat_consensus_id, repeat_consensus " .
			  "from repeat_consensus where repeat_class = \"trf\" ".
			  "AND length(repeat_consensus) = $i" );

  $sth->execute();

  my ( %rc_2_norm, %rc_2_seq, %norm_2_rc );

  while( my $arr = $sth->fetchrow_arrayref() ) {
    my ( $rc_id, $rc_seq ) = @$arr;
    $all_cons_count++;
    my $norm_seq = norm_seq( $rc_seq );
    $rc_2_norm{ $rc_id } = $norm_seq;
    $rc_2_seq{ $rc_id } = $rc_seq;
    if( $norm_seq eq $rc_seq ) {
      $norm_2_rc{ $norm_seq } = $rc_id;
    }
  }

  # iterate through the rc_ids and write old tab new lines out
  for my $rc_id ( keys %rc_2_norm ) {
    my $norm = $rc_2_norm{ $rc_id };
    if( exists $norm_2_rc{ $norm } ) {
      my $new_rc = $norm_2_rc{ $norm };
      if( $rc_id != $new_rc ) {
	print FH "$rc_id\t$new_rc\n";
      }
    } else {
      $norm_2_rc{ $norm } = $rc_id;
    }
  }
}

my $sth = $db->prepare( "SELECT repeat_consensus_id, length( repeat_consensus ) " .
			"FROM repeat_consensus WHERE repeat_class = \"trf\" " .
			"AND length(repeat_consensus) > 8" );

$sth->execute();

my %len_2_rc;
my $length_removal_count = 0;

while( my $arr = $sth->fetchrow_arrayref() ) {
  my ( $rc_id, $length ) = @$arr;
  if( exists $len_2_rc{$length} ) {
    # map this rc to the first one of that length
    print FH "$rc_id\t".$len_2_rc{ $length }."\n";

  } else {
    $len_2_rc{ $length } = $rc_id;
  }
}

print STDERR "File written ";

close( FH );
$db->do( "create table tmp_old_new_rcid ( old_id int not null, new_id int not null, key old_idx( old_id))" );

$db->do( "load data local infile '$filename' into table tmp_old_new_rcid" );

print STDERR "and uploaded.\n";

$db->do( "update repeat_feature rf, tmp_old_new_rcid tonr " .
	 "set rf.repeat_consensus_id = tonr.new_id " .
	 "where rf.repeat_consensus_id = tonr.old_id " );

print STDERR "Repeat_features updated.\n";

$db->do( "delete from rc " .
	 "using repeat_consensus rc, tmp_old_new_rcid tonr " .
	 "where rc.repeat_consensus_id = tonr.old_id " );

$db->do( "update repeat_consensus " .
	 "set repeat_consensus = concat( length( repeat_consensus ), \"(N)\") " .
	 "where repeat_class = \"trf\" " .
	 "and length( repeat_consensus ) > 8" );

$db->do( "drop table tmp_old_new_rcid" );

unlink( $filename );

exit;



# find a representation of the input sequence that is canonical
# rotate and revcomp all possibilities and take alphabetical lowest
sub norm_seq {
  my $seq = shift;
  my $test_seq = $seq;
  
  my @all_seq = ();
  push( @all_seq, $seq );

  for( my $i=1; $i<2*length( $seq ); $i++ ) {
    if( $i != length( $seq )) {
      $test_seq = rotate( $test_seq );
    } else {
      $test_seq = rev_comp( $test_seq );
    }
    my $new_seq = $test_seq;
    push( @all_seq, $new_seq );
  }

#  print "old: ",join( " ", @all_seq),"\n";
  @all_seq = sort {$a cmp $b} @all_seq;
#  print "new: ",join( " ", @all_seq),"\n";
  return $all_seq[0];
}
  
	       

sub rotate {
  my $string = shift;
  return substr( $string, 1 ) . substr( $string, 0, 1 );
}

sub rev_comp {
  my $string = shift;
  $string = reverse( $string );
  $string =~
    tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  return $string;
}

sub usage {

  print STDERR <<EOC

  Usage: perl shrink_trfs.pl -user .. -port .. -pass .. -host .. -dbname ..
         Tries to minimize the number of repeat_consensi by only having one 
         trf consensus for each length of repeat_consensus longer than 8. 
         Shorter repeat consensi are reduced to the alphabetical minimum
         of all rotate / revcomp equivalent repeats.
  
EOC
;
  exit;
}

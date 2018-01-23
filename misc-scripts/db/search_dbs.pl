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


# this script will lokk through all listed databases for 
# database names which match given regexp.
# This version should be ok if one or more of the given database
# servers is not available.


use strict;
use warnings;

use DBI;

# seconds before we give up on database
my $timeout = 5;

my @all_locators = qw {
		     genebuild1:3306:ensro::
		     genebuild2:3306:ensro::
		     genebuild3:3306:ensro::
		     genebuild4:3306:ensro::
		     genebuild5:3306:ensro::
		     genebuild6:3306:ensro::
		     genebuild6:3307:ensro::

		     compara1:3306:ensro::
		     compara2:3306:ensro::
		     compara3:3306:ensro::

		     mart1:3306:ensro::
		     mart2:3306:ensro::

		     ens-genomics1:3306:ensro::
		     ens-genomics2:3306:ensro::

		     ens-research:3306:ensro::
		     ens-research:3309:ensro::

                     ensdb-1-11:3319:ensro::
                     ensdb-1-11:3317:ensro::

		     ens-staging:3306:ensro::
		     ens-livemirror:3306:ensro::

		     };


my $pattern=shift;
my $size = 0;
my %sizes;

if( @ARGV ) {
  $size = 1;
}

if( !$pattern ) {
  printf( "You need to supply a regexp as argument.\n" );
  printf( "Use -list if you want to see which databases are configured.\n" );
  printf( "If you have any parameters after the regexp, the program will print\n" );
  printf( " database sizes as well.\n" );
  printf( "\nExample: search_dbs.pl \"^stabenau\" -sizes | sort -knr3\n" ); 
  exit;
}

my $list = 0;

if( $pattern eq "-list" ) {
  $list = 1;
}


for my $loc ( @all_locators ) {
  my @elems = split( ":", $loc );
  my @dbnames = ();
  %sizes = ();
  my $dsn = sprintf( "DBI:mysql:host=%s;port=%s", $elems[0], $elems[1] );

  $SIG{ALRM} = sub{die("timeout");};

  eval {
    alarm $timeout;
    my $db = DBI->connect( $dsn, $elems[2], $elems[3], { RaiseError => 1 } );
    my $res = $db->selectall_arrayref( "show databases" );
    @dbnames = map { $_->[0] } @$res;
    $db->disconnect();
    alarm 0;
  };

  if( !@dbnames ) {
    print STDERR "$loc NOT OK\n";
  } else {
    if( $size ) {
      for my $dbname ( @dbnames ) {
	if( $dbname =~ /$pattern/ ) {
	  eval {
	    alarm $timeout;
	    my $db = DBI->connect( $dsn, $elems[2], $elems[3], { RaiseError => 1, PrintError=> 0 } );
	    $db->do( "use $dbname" );
	    my $t_status = $db->selectall_arrayref( "show table status" );
	    my $size = 0;
	    map { $size += $_->[6]; $size += $_->[8] } @$t_status;
	    print "$loc $dbname $size\n";
	    $db->disconnect();
	    alarm 0;
	  };
	  if( $@ ) {
	    print( "Problem on $loc $dbname.\n ", $@ );
	  }
	}
      }
    } else {
      if( $list ) {
	print STDERR "$loc ok\n";
      } else {
	for my $dbname ( @dbnames ) {
	  if( $dbname =~ /$pattern/ ) {
	    print "$loc $dbname\n";
	  }
	}
      }
    }
  }
}

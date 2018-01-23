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


###############################################################################################################
# multidb_sql script
# given a pattern of databases to query, runs a query accross all databases
# can take a file as input, to apply sql without printing out results
# USAGE : perl multidb_sql.pl -dbpattern _core_ -expr "select count(*) from gene"
################################################################################################################




use strict;
use warnings;

use DBI;
use Getopt::Long;

my ( @hosts, $host, $user, $pass, $port, $expression, $dbpattern, $file, $result_only ); 

$user = 'ensro' ;  
$pass = '' ;
$port = 4519;
@hosts = qw ( mysql-ens-sta-1 ) ;

GetOptions( "host|dbhost=s", \$host,
            "hosts|dbhosts=s", \@hosts,
	    "user|dbuser=s", \$user,
	    "pass|dbpass=s", \$pass,
	    "port|dbport=i", \$port,
	    "expr=s", \$expression,
	    "file=s", \$file,
	    "dbpattern|pattern=s", \$dbpattern,
            "result_only!", \$result_only,
	  );

if ($host) {
  @hosts = $host ;
}
else {
  @hosts = split(/,/,join(',',@hosts)) ;
}

foreach my $host ( @hosts ) {
  my $dsn = "DBI:mysql:host=$host";
  if( $port ) {
    $dsn .= ";port=$port";
  }

  my @expressions;

  if( $file ) {
    local *FH;
    if( ! -r $file ) {
      die ( "File $file not readable" );
    }
    open( FH, "<$file" );
    my $exp;
    while( my $line = <FH> ) {
      if( $line =~ /;$/ ) {
        $line =~ s/;$//;
        $exp .= " ".$line;
        push( @expressions, $exp );
        $exp = "";
      } else {
        $exp .= " ".$line;
      }
    }
    if( $exp ) {
      push( @expressions, $exp );
    }  
    close FH;
  }

  my $db = DBI->connect( $dsn, $user, $pass ) or
      die "Unable to connect to database(s): $dsn, $user";
  my @dbnames = map {$_->[0] } @{ $db->selectall_arrayref( "show databases" ) };
  for my $dbname ( @dbnames) {   
    next if ( $dbname !~ /$dbpattern/ );
    unless ($result_only) { 
      print  "$dbname" . "@" . "$host\n";
    }
    if(( ! $expression ) && ( !$file )) {
      next;
    }
  
    $db->do( "use $dbname" );
    if( $file ) {
      for my $sql ( @expressions ) {
        print STDERR "Do $sql\n";
        $db->do( $sql );
      }  
    } elsif( $expression =~ /^\s*select/i || $expression =~ /^\s*show/i || $expression =~ /^\s*desc/i ) {
      my $res = $db->selectall_arrayref( $expression );
      my @results = map { join( " ", @$_ ) } @$res ;
      my $db_name_off = 0 ;
      for my $result ( @results ) {
        if($result_only){ 
          unless ($db_name_off){
            $db_name_off =1 ;
            print  "==> $dbname" . "@" . $host . " :\n";
          }
        }
        print  "    Result: ",$result,"\n";
      }
    } else {
      $db->do( $expression );
      print  "  done.\n";
    }
  }  
}

sub usage {
  print STDERR <<EOF

             Usage: multidb_sql options
 Where options are: -host hostname 
                    -user username 
                    -pass password 
                    -port port_of_server optional
                    -dbpattern regular expression that the database name has to match
                    -expr sql statement you want to execute. 
                          if omitted, just print database names matching
                          if select, show or describe prints results
                    -file Apply sql in file to all databases. Does not print results.
                    -result_only only prints out databases for which results where found

EOF
;
  exit;
}


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
# convert_tables_MyISAM_InnoDB.pl script
# given a pattern of databases to query, lists all tables and whether they are InnoDB or MyISAM
# can convert a list of tables, or all tables, to MyISAM, or to InnoDB, as specified
# USAGE : perl convert_tables_MyISAM_InnoDB.pl -dbpattern _core_ -convert_to MyISAM -convert_all
################################################################################################################

# $Source: /cvsroot/ensembl/ensembl-personal/genebuilders/scripts/convert_tables_MyISAM_InnoDB.pl,v $
# $Revision: 1.4 $


use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use DBI; 

my ($host, $user, $pass, $port, $dbpattern, $verbose, $convert_all, $convert_to, $tables);

$user = 'ensro';
$pass = '';
$port = 3306;
$host = 'ens-staging';
$verbose = 0;
$dbpattern = '';

GetOptions( "host|h|dbhost=s"     => \$host,
            "user|u|dbuser=s"     => \$user,
            "pass|p|dbpass=s"     => \$pass,
            "port|P|dbport=i"     => \$port,
            "dbpattern|pattern=s" => \$dbpattern,
            'verbose!'            => \$verbose,
            'convert_all!'        => \$convert_all,  # flag to convert all tables
            'convert_to=s'        => \$convert_to,   # engine to covert the tables to ( MyISAM, InnoDB )
            'tables=s'            => \$tables,       # comma-separated list of tables to convert ie exon,exon_transcript
          );

my $dsn = "DBI:mysql:host=$host";
if( $port ) {
  $dsn .= ";port=$port";
}
my $db = DBI->connect( $dsn, $user, $pass) ;    
my @dbnames = map {$_->[0] } @{ $db->selectall_arrayref( "show databases" ) };
for my $dbname ( @dbnames ) {
  if( $dbname =~ /$dbpattern/ ) {     
    print "connecting to $dbname\n" ; 

    my $dsn_info = sprintf( "DBI:mysql:database=%s;host=%s;port=%s", 'information_schema', $host, $port)  ;  
    my $db_info = DBI->connect( $dsn_info, $user, $pass ) ;   

    my %engine_tables = %{  check_for_tables($db_info, $dbname, $verbose) } ;   
    my @tables_to_convert ;   
      
    if ( $tables ) { 
      @tables_to_convert = split /\,/, $tables ;   
    } elsif ( $convert_all ) { 
      if ( $convert_to =~m/MyISAM/) {    
        die('There is no table to convert! Exiting...'."\n") unless (exists $engine_tables{"InnoDB"}); 
        @tables_to_convert = @{$engine_tables{"InnoDB"}};
      }elsif ( $convert_to =~ m/InnoDB/ ) { 
        die('There is no table to convert! Exiting...'."\n") unless (exists $engine_tables{"MyISAM"}); 
        @tables_to_convert = @{$engine_tables{"MyISAM"}} ;  
      }   
    }
   
    unless ( $convert_to ) {  
      print "\n\n\nTo convert selected tables to use a different storage engines, use : \n\n"     ; 
      print "\t\t-convert_to [MyISAM|InnoDB]' -tables exon_transcript,job,exon\n\n\n"
           ." OR convert all tables with :\n\n\t\t-convert_to MyISAM -convert_all  \n" ;   
    } else {  

      print "\n\nWill convert these tables : \n" ;  
      for my $table  ( @tables_to_convert ) { 
          printf "%-10s  ===>  %-10s\n", $table, $convert_to ;  
      }  
      print "\n\n\tARE YOU SURE YOU WANT TO GO AHEAD ??? ( Y / N )" ; 
      if ( get_input_arg() ) {
        convert_all_tables ( $dbname, \@tables_to_convert, $convert_to   ) ;  
      } else {  
        print "Stopping as you don't want to go ahead. no tables have been converted \n" ; 
      } 
     } 
    $db_info->disconnect() ; 
  }
}  


sub convert_all_tables {   
  my ( $dbname, $tables_to_convert, $convert_to )  = @_ ;   

   my $dsn_info = sprintf( "DBI:mysql:database=%s;host=%s;port=%s", $dbname, $host, $port)  ;  
   my $dbh  = DBI->connect( $dsn_info, $user, $pass ) ;      

   for my $table ( @$tables_to_convert ) { 
        my $sql= "alter table $table engine =\'".$convert_to."\';"  ; 
        $dbh->do($sql) ;  
        print "table $table converted to $convert_to\n" ; 
   }      
   $dbh->disconnect(); 
} 
 
sub check_for_tables { 
   my ($db_info, $dbname, $verbose) = @_;   

   my $sql ="select table_type, table_name, engine from tables where table_schema =\'".$dbname. "\' " ; 

   my $lines = $db_info->selectall_arrayref($sql) ;   

   my %engine_types; 

   for my $r( @$lines) {   
     my @rows  = @$r ; 
     if ($rows[0] eq 'BASE TABLE') {
       push @{ $engine_types{$rows[2]}}, $rows[1] ; 
     } elsif ($rows[0] eq 'VIEW') {
       push @{ $engine_types{$rows[0]}}, $rows[1] ;
     }
   }  
   for ( keys %engine_types )  {  
     print uc($_) . "  " . scalar(@{$engine_types{$_}}) . " tables $dbname\n" ;   
   }  

   if ( $verbose ) {     
     print "\n\nTable types found :\n---------------------------------------\n\n" ; 
     for my $engine( keys %engine_types )  {    
        for my $table ( @{ $engine_types{$engine}} ) {  
          printf "%-10s%-10s\n", $engine, $table ; 
        }  
     }  
   } 
   return \%engine_types ; 
}  


=head2 get_input_arg  ( Bio::EnsEMBL::Analysis::Tools::Utilities) 

  Function  : waits for input from STDIN and returns '1' if input =~m/y/i         
              and '0' if input matches /n/i.  
  Returntype: 1 or 0 
  Exceptions: none

=cut

sub get_input_arg {
  while (defined (my $line=<STDIN>)){
   chomp($line) ;
   if ( $line=~m/y/i){
      return 1 ;
   }elsif( $line =~m/n/i){
     return 0 ;
   }
   print "Wrong input - only answer 'y' or 'n'\n" ;
  }
}



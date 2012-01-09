#!/usr/bin/env perl
##########################################
#
# Simple but handy script that generate the input file
# neded by the script CopyDBoverServer.pl
#
# Benoit Ballester <benoit@ebi.ac.uk>
#
# Please post comments/questions to the Ensembl development list
# <dev@ensembl.org>
#
#########################################
use strict;
use warnings;
use DBI;
use Getopt::Long;

my $help= 0;
my $sourceHost="ens-staging";
my $sourcePort="3306";
my $sourceUser="";
my $sourcePwd="";
my $destinationHost="mart2";
my $destinationPort="3306";
my $limit='';
my $target_location = '';

my $usage = "\nUsage: $0 -sourceHost mart1 -sourceUser xxx -destinationHost mart2  -limit %42%\n
  -help or -h      [for help]

  -sourceHost      [default: ens-staging]
  -sourcePort      [default: 3306       ]
  -sourceUser      [default:            ]
  -destinationHost [default: mart2      ]
  -destinationPort [default: 3306       ]
  -limit           [eg. %core_42%       ]
  -target_location [default: blank if you want standard data locations]\n

The limit option will limit the databases being copied according to your limit criteria. 
With   -limit %core_42%   only ensembl core 42 databases will be copied\n\n
";

GetOptions('help|h' => \$help,
			 'sourceHost=s' => \$sourceHost,
			 'sourcePort=s' => \$sourcePort,
			 'sourceUser=s' => \$sourceUser,
			 'sourcePwd=s' => \$sourcePwd,
			 'destinationHost=s' => \$destinationHost,
			 'destinationPort=s' => \$destinationPort,
			 'target_location=s' => \$target_location,
			 'limit=s' => \$limit);

if ($help || scalar @ARGV == 0 ) {
#if ($help) {
    print $usage;
  exit 0;
}
#--------------connect to MySQL Source Host
my $dbh   = DBI->connect ("DBI:mysql:host=$sourceHost:port=$sourcePort",
                          $sourceUser,
                          $sourcePwd)
    or die "Can\'t connect to database: $DBI::errstr\n";

#--------------prepare and execute the query
my $sql;
if (!$limit){
    $sql = "show databases ;";
}else {
    $sql = "show databases like \"".$limit."\" ;";
}

my $sth = $dbh->prepare( $sql);

$sth->execute( );
while ( my @row = $sth->fetchrow ){
    
    my $result = sprintf ("%s %d %50s %s %d %s",$sourceHost.".internal.sanger.ac.uk", $sourcePort, $row[0], $destinationHost.".internal.sanger.ac.uk ", $destinationPort, $row[0]);
    $result .= "  $target_location" if $target_location;
    print $result . "\n";
    
}

$sth->finish( );

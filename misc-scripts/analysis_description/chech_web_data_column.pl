#!/usr/local/ensembl/bin/perl

=pod

=head1 NAME

  check_web_data_column.pl

=head1 SYNOPSIS

  script to check that the web_data column can be eval'd into a hash.

=head1 DESCRIPTION

   This column has been giving problems due to the peculiarity of the data 
   (string containing quotes that has to be eval'd into a hash ref). In 
   order to check it, this script should be run after load_analysis_description
   and apply_rules

=head1 OPTIONS

  Database options

    -dbhost      host name for database (default=ens-staging)
    -dbport      For RDBs, what port to connect to (default=3306)
    -dbname      For RDBs, what name to connect to (dbname= in locator)
    -pattern     check databases matching this PATTERN
                 Note that this is a database pattern of the form %core_53_%
                 rather than a regular expression
    -dbuser      For RDBs, what username to connect as (dbuser= in locator)
    -dbpass      For RDBs, what password to use (dbpass= in locator)
    -help print out documentation

=head1 EXAMPLES

  In order to run it for a set of databases

  perl check_web_data_column -dbuser ensro -pattern '%core_53%_'

=cut

use strict;
use warnings;
use Getopt::Long;
use DBI;
use DBD::mysql;
use Data::Dumper;

my ($dsn,$dbh);

my $dbhost = 'ens-staging';
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;
my $pattern;
my $help = 0;

&GetOptions(
	    'host|dbhost=s'      => \$dbhost,
	    'port|dbport=s'      => \$dbport,
	    'user|dbuser=s'      => \$dbuser,
	    'pass|dbpass=s'      => \$dbpass,
            'dbname=s'           => \$dbname,
            'pattern=s'          => \$pattern,
            'help|h!'            => \$help
);

if (!$dbhost){
    print ("Need to pass a dbhost\n");
    $help =1;
}
if (!$dbname and !$pattern){
    print("Need to enter either a database name in -dbname or a pattern in -pattern");
    $help = 1;
}

if ($help){
    usage();
}

#connect to database
$dsn = "DBI:mysql:host=" . $dbhost . ";port=" . $dbport;

eval{ 
  $dbh = DBI->connect($dsn, $dbuser, $dbpass, 
		      {'RaiseError' => 1,
		       'PrintError' => 0});
};

# get all database names that match pattern
my ($sth, $sql);
my $sql_pattern = $pattern || $dbname;
$sql = "SHOW DATABASES LIKE '". $sql_pattern ."'";
$sth = $dbh->prepare($sql);
$sth->execute;

my $analysis_id;
my $display_label;
my $web_data;
my $ref_web_data;
my $db_found = 0;
#for each of the database, check the web_data in the analysis_description
while (my ($dbname) = $sth->fetchrow_array) {
    $sql = qq{SELECT analysis_id, display_label, web_data from $dbname.analysis_description};
    my $sth1 = $dbh->prepare($sql);
    $sth1->execute;
    $sth1->bind_col(1,\$analysis_id);
    $sth1->bind_col(2,\$display_label);
    $sth1->bind_col(3,\$web_data);

    while ($sth1->fetch){
	if (defined $web_data and $web_data ne ''){
	    $ref_web_data = eval($web_data); 
	    #print Dumper($ref_web_data);
	    if (ref($ref_web_data) ne 'HASH'){
		print "Analysis $analysis_id with display_label $display_label in $dbname has a wrong web_data column--$web_data-- cannot be eval into hash\n";
		$db_found = 1;
	    }
	}
    }
}
print "Web data column in $sql_pattern looks right\n" if (!$db_found);

sub usage{
    exec('perldoc',$0);
    exit;
}

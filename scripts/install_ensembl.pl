#!/usr/local/bin/perl

use strict;

BEGIN {
    unshift(@INC,"../modules");
}


#
# Ok. We have to be super paranoid and very nice
# to our users here.
# 


#
# say hello
#

print STDERR "Welcome to Ensembl!\n\n";

#
# See if we can load DBI 
#

eval {
    require DBI;
};

if( $@ ) {
    print STDERR "I am sorry, I was unable to find the DBI (database interface) for Perl\n";
    print STDERR "You need to install a SQL database and the DBI interface, for example, MySQL\n";
    print STDERR "and MySQL DBI. The DBI interfaces are kept on CPAN (not the MySQL site)\n";
    exit(0);
}

print STDERR "... DBI checks out\n";

#
# See if we can load the bioperl root object and AnnSeq object
# AnnSeq is 0.6 or later...
#

eval {
    require Bio::Root::Object;
    require Bio::AnnSeq;
};

if( $@ ) {
    print STDERR "I am sorry, I was unable to find the bioperl root object\n";
    print STDERR "You need to have a bioperl 0.6 release installed\n";
    print STDERR "You can one from ftp://ftp.sanger.ac.uk/pub/ensembl/software\n";
    exit(0);
}
print STDERR "... bioperl checks out\n";

#
# See if we can load ensembl root object
#

eval {
    require Bio::EnsEMBL::DBSQL::Obj;
};

if( $@ ) {
    print STDERR "I am sorry, I was unable to load the Ensembl root object\n";
    print STDERR "You might be running this somewhere other than the ensembl software\n";
    print STDERR "distribution, in which case, please make sure your Perl include path has\n";
    print STDERR "ensembl modules in it\n\n";
    print STDERR "Another possibility is that you have an old version of bioperl\n";
    print STDERR "You need to have a bioperl 0.6 release installed\n";
    print STDERR "You can one from ftp://ftp.sanger.ac.uk/pub/ensembl/software\n";
    print STDERR "\n\nActual exception $@\n";
    exit(0);
}

print STDERR "... ensembl checks out\n";

#
# Tell our user that he/she is looking good
#

print STDERR "The ensembl software looks ready to run!\n\n";
print STDERR "Looking at your DBI installation. I'll need some details about this to establish a connection\n\n";

print STDERR "What is the name of the ensembl database [ensdev]\n\n(If you haven't created a database, no worries\nSimply Cntr-C out of this and go\nmysqladmin -u username create ensdev\nYou may need to talk to whoever install mysql to get permissions to create ensdev)\n\nWhat is the name of the ensembl database [ensdev]";

my $db = <STDIN>;
chomp $db;
if( $db =~ /^\s*$/ ) {
    $db = 'ensdev';
}

#
# Ask whether it is ok to use a MySQL driver?
#

print STDERR "What username shall I try to connect to the MySQL db [root]?\n";
my $user = <STDIN>;
chomp $user;
if( $user =~ /^\s*$/ ) {
    $user = 'root';
}

print STDERR "What password shall I try to connect to the MySQL db []?\n";
my $pass = <STDIN>;
chomp $pass;
if( $pass =~ /^\s*$/ ) {
    $pass = '';
}

print STDERR "What host shall I try to connect to the MySQL db [localhost]?\n";
my $host = <STDIN>;
chomp $host;
if( $host =~ /^\s*$/ ) {
    $host = 'localhost';
}

my $dsn;
my $dbh;
eval {
    $dsn = "DBI:mysql:database=$db;host=$host;";

    $dbh = DBI->connect("$dsn","$user",$pass);

};

if( $@ || ! defined $dbh ) {
    print STDERR "Well. That did not work. I cannot connect with DBI to the mysql database\n";
    print STDERR "I was using the following connection string [$dsn]\nThe exception thrown was\n$@\n\n";
    print STDERR "I'm afraid you need to debug this before continuing\n";
    exit(0);
}

$dbh->disconnect();
$dbh = 0;

#
# Find the DNA and table data 
#

my $table_file;

if( -e "../data/ensembl.sql" ) {
    $table_file = "../data/ensembl.sql";
}

print STDERR "Can you locate the ensembl.sql file [$table_file]?";

do {
    my $fname = <STDIN>;
    chomp $fname;
    if( $fname =~ /^\s*$/ && defined $table_file ) {
	last;
    }

    if(! -e $fname ) {
	print STDERR "Cannot open $fname\n";
    } else {
	$table_file = $fname;
    }

} while ( ! defined $table_file );

my $dna_file;



#
# Sanity checks
#

if ( !-e "../sql/table.sql" || !-e $table_file ) {
    print STDERR "Cannot find one of table.seq, dna file or table data file. Sorry.\n\n";
    exit(0);
}


#
# Ok. Reinitialise the db?
#

print STDERR "I am now going to intialize the database.\n";
print STDERR "and then reload it with new data\n\n";
print STDERR "This will trash any current data in ensdev\n";
print STDERR "Please answer yes or no to this reintialisation\n";

my $q = <STDIN>;
chomp $q;

if( $q ne 'yes' ) {
    print STDERR "\n\nBye...\n";
}

print "... reinitialisation of tables\n";

system("mysql -h $host -u $user $db < ../sql/table.sql");

print ".......done\n";

print "... loading data\n";

system("mysql -h $host -u $user $db < $table_file");

print ".......done\n";

print STDERR "Your distribution is read to run!\n";






#!/usr/local/bin/perl

=head1 NAME - gtf dump

    Dumps EnsEMBL genes in GTF format (Gene Transfer Format)

=head1 SYNOPSIS - 

    gtfdump -getall -dbname test -dumpfile all_genes.gtf

    gtfdump -getall -clone_file new_clones.list -dumpfile new_genes.gtf

    gtfdump -dumpfile AB019437-genes.gtf AB019437

=head1 DESCRIPTION

    This module dumps EnsEMBL genes from a database to a file in GTF
    format. The dump happens on a clone by clone basis. The clones to be
    dumped can be chosen either by specifying the -getall option, which dumps
    all genes from a database, or by using the -clonelist and providing a list    of clones to dump, or by simply passing the clone ids on the command lines    after all the other options.

    The actual dumping happens in the Bio::EnsEMBL::Utils::GTF_handler module,    which also handles the parsing of GTF files.

=head1 OPTIONS

    -dbtype   database type (only needed for TimDB)

    -host     host name for the database (gets put as host= in locator) 

    -port     for RDBs, what port to connect to for the "to" database (port= in locator)

    -dbname   for RDBs, what database name to connect to (dbname= in locator)

    -dbuser   for RDBs, what username to connect as to the database (dbuser= in locator)

    -dbpass   for RDBs, what password to use to connect to to the database (dbpass= in locator)

    -module   module name to load to (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -getall   all clones from the database [not applicable to timdb]

    -clonelist  read in on stdin a list of clones, one clone per line

    -dumpfile  name of the file to dump genes to 
    
    -help      Displays this documentation with PERLDOC

=cut

use Bio::EnsEMBL::Utils::GTF_handler;
use Bio::EnsEMBL::DBLoader;
use strict;
use Getopt::Long;

#Database options
my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'f15';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';

#Clone options
my $getall;
my $clonelist;

#Other options
my $dumpfile;
my $help;

&GetOptions( 
	     'dbtype:s'  => \$dbtype,
	     'host:s'    => \$host,
	     'port:n'    => \$port,
	     'dbname:s'  => \$dbname, 
	     'dbuser:s'  => \$dbuser,
	     'dbpass:s'  => \$dbpass,
	     'module:s'  => \$module,
             'getall'    => \$getall,
             'clonelist' => \$clonelist,
	     'dumpfile:s'=> \$dumpfile,
	     'h|help'    => \$help
	     );

my $gtfh=Bio::EnsEMBL::Utils::GTF_handler->new();

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
print STDERR "Using $locator for todb\n";

my $db =  Bio::EnsEMBL::DBLoader->new($locator);
my @clones;

if ($clonelist) {
   while( <> ) {
   my ($en) = split;
   push(@clones,$en);
   }
}
elsif ( $getall) {
   @clones = $db->get_all_Clone_id();
   print STDERR scalar(@clones)." clones found in DB\n";
}
else {
   @clones = @ARGV;
}
if ( $getall) {
   @clones = $db->get_all_Clone_id();
   print STDERR scalar(@clones)." clones found in DB\n";
}

#push @clones,'AB019437';

foreach my $clone_id (@clones) {
    my $clone=$db->get_Clone($clone_id);
    my @genes=$clone->get_all_Genes();
    $gtfh->dump_genes($dumpfile,@genes);
}




#!/usr/local/bin/perl

=head1 NAME - gtf dump

    Dumps EnsEMBL genes in GTF format (Gene Transfer Format)

=head1 SYNOPSIS - 

    gtfdump [ < ids-file | -getall]  [ -clones | -chrs]  -dbname test [ -dumpfile all_genes.gtf | -separefiles ]  [ ids ... ]

=head1 DESCRIPTION

    This module dumps EnsEMBL genes from a database to a file in GTF
    format. The dump happens on a clone by clone or chr by chr basis. The
    things to be dumped can be chosen either by specifying the -getall
    option, which dumps genes  on all things in the database, or 
    reading the list of ids from stdin, or by simply
    passing the clone ids on the command lines after all the other
    options.

    The actual dumping happens in the Bio::EnsEMBL::Utils::GTF_handler
    module, which also handles the parsing of GTF files.

=head1 OPTIONS

    -dbtype   database type (only needed for TimDB)

    -host     host name for the database (gets put as host= in locator) 

    -port     for RDBs, what port to connect to for the "to" database (port= in locator)

    -dbname   for RDBs, what database name to connect to (dbname= in locator)

    -dbuser   for RDBs, what username to connect as to the database (dbuser= in locator)

    -dbpass   for RDBs, what password to use to connect to to the database (dbpass= in locator)

    -module   module name to load to (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -getall   all clones or chrs from the database 

    -separatefiles genes of each clone or chr go to separate <the_id>.gtf

    -dumpfile  name of the file to dump all genes to  (stdout if not specified)
    
    
    -help      Displays this documentation with PERLDOC

=cut

use Bio::EnsEMBL::Utils::GTF_handler;
use Bio::EnsEMBL::DBLoader;
use strict;
use Getopt::Long;

#Database options
my $dbtype = 'rdb';
my $host   = 'kaka.sanger.ac.uk';
my $port   = undef; # '410000';
my $dbname = 'ensembl';
my $dbuser = 'anonymous';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';

#Clone options
my $getall;
my $clones;
my $chrs;

#Other options
my $dumpfile;
my $separatefiles;
my $help;

&GetOptions( 
	     'dbtype:s'  => \$dbtype,
	     'host:s'    => \$host,
	     'port:n'    => \$port,
	     'dbname:s'  => \$dbname, 
	     'dbuser:s'  => \$dbuser,
	     'dbpass:s'  => \$dbpass,
	     'module:s'  => \$module,
             'clones'    => \$clones,
             'chrs'    => \$chrs,
             'getall'    => \$getall,
	     'dumpfile:s'=> \$dumpfile,
             'separatefiles' =>\$separatefiles,
	     'h|help'    => \$help
	     );
my $usage = "gtfdump [ < ids-file | -getall]  [ -clones | -chrs]  -dbname test [ -dumpfile all_genes.gtf | -separefiles ]  [ ids ... ]";

die $usage,"\n" if  $help;

if ( !($clones || $chrs) || ($clones && $chrs )) {
    die "need to specify exactly one of [ -clones | -chrs ]\n";
}

if ($dumpfile && $separatefiles) {
    die "specify at most one of [ -dumpfile FILE | -separatefiles ]\n";
}

if ( $dumpfile && !$separatefiles) {
    open (DUMP,">$dumpfile") || die ("Could not dump to $dumpfile\n");
} else {
    *DUMP=*STDOUT;
}

my $gtfh=Bio::EnsEMBL::Utils::GTF_handler->new();

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
print STDERR "Using $locator for todb\n";

my $db =  Bio::EnsEMBL::DBLoader->new($locator);
my @ids;

if ( $getall) {
    if ($clones) {
        @ids = $db->get_all_Clone_id();
    } else {
        @ids = $db->get_all_chr_ids('UCSC');  
    }
    warn scalar(@ids)." ids found in DB\n";
} elsif ( @ARGV > 0 ) {
    @ids = @ARGV;
} else {                                # lines on stdin
    while( <> ) {
        my ($en) = split;
        push(@ids,$en);
    }
}

my $n;
my $stap = $db->get_StaticGoldenPathAdaptor();

foreach my $id (@ids) {
    warn "doing $id ...\n";
    my @genes;
    if ($clones) {
        my $clone=$db->get_Clone($id);
        @genes=$clone->get_all_Genes();
    } else {
        my $vc = $stap->fetch_VirtualContig_by_chr_name($id);
        @genes = $vc->get_all_Genes();
    }
        

    next if int(@genes) == 0;
    warn scalar(@genes)." genes found in db\n";

    if ($separatefiles) {
        my $file = "$id.gtf";
        open(F,">$id.gtf") || die "Cannot open $file: $!";
        $gtfh->dump_genes(\*F,@genes);
        close(F);
    } else { 
        $gtfh->dump_genes(\*DUMP,@genes);
    }

#     last if $n++ > 100;
}




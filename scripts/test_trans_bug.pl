#!/usr/local/bin/perl

=head1 NAME

trans2xml

=head1 SYNOPSIS
 
 trans2xml.pl

=head1 DESCRIPTION

This script produces GAME XML file outputs of the whole database, where each contig produces a new
annotation ID and each transcript being a new feature_set annotation.

=head1 OPTIONS

    -dbtype    Database tpye (only used for TimDB)

    -host      host name for database (gets put as host= in locator)

    -port      For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -help      Displays script documentation with PERLDOC

=cut

use strict;
use Bio::EnsEMBL::DBLoader;
use Bio::SeqIO;
use Getopt::Long;

my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'test_trans';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $getall;
my $usefile;
my $help;

&GetOptions( 	     
	     'host:s'   => \$host,
	     'port:n'   => \$port,
	     'dbname:s' => \$dbname,
	     'dbuser:s' => \$dbuser,
	     'dbpass:s' => \$dbpass,
	     'module:s' => \$module,
	     'h|help'   => \$help,
	     'getall'   => \$getall,
	     'usefile'=> \$usefile
	     );

if ($help) {
    exec('perldoc', $0);
}

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);
my $seqio;
my @trans;

if ($help){
    exec('perldoc', $0);
}

if( $usefile == 1 ) {
    while( <> ) {
	my ($en) = split;
	push(@trans,$en);
    }
}
else {
    @trans = @ARGV;
}

foreach my $trans_id (@trans) {
    print STDERR "\nDumping transcript $trans_id\n";
    my $gene_obj=Bio::EnsEMBL::DBSQL::Gene_Obj->new($db);
    my $trans = $gene_obj->get_Transcript($trans_id);
    print STDERR "First Exon: ".$trans->first_exon->id."\n";
    print STDERR "Last Exon: ".$trans->last_exon->id."\n\n";
    print STDERR "This trasncript contains the following exons:\n";
    foreach my $exon ($trans->each_Exon) {
	print STDERR "      id: ".$exon->id."\n";
	print STDERR "  strand: ".$exon->strand."\n";
	print STDERR "   start: ".$exon->start."\n";
	print STDERR "     end: ".$exon->end."\n";
    }

    print STDERR "Now redumping in vc coordinates\n";
    my $trans=$db->get_Transcript_in_VC_coordinates($trans_id);
    print STDERR "\nDumping transcript ".$trans->id."\n";
    print STDERR "First Exon: ".$trans->first_exon->id."\n";
    print STDERR "Last Exon: ".$trans->last_exon->id."\n\n";
    print STDERR "This trasncript contains the following exons:\n";
    foreach my $exon ($trans->each_Exon) {
	print STDERR "      id: ".$exon->id."\n";
	print STDERR "  strand: ".$exon->strand."\n";
	print STDERR "   start: ".$exon->start."\n";
	print STDERR "     end: ".$exon->end."\n";
    }
}

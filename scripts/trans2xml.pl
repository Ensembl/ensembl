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
use Bio::EnsEMBL::TimDB::Obj;
use Bio::SeqIO;
use Getopt::Long;

my $host   = 'sol28';
my $port   = '410000';
my $dbname = 'ensdev';
my $dbuser = 'ensembl';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBOLD::Obj';
my $help;

&GetOptions( 	     
	     'host:s'   => \$host,
	     'port:n'   => \$port,
	     'dbname:s' => \$dbname,
	     'dbuser:s' => \$dbuser,
	     'dbpass:s' => \$dbpass,
	     'module:s' => \$module,
	     'h|help'   => \$help,
	     );

if ($help) {
    exec('perldoc', $0);
}

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);
my $seqio;

foreach my $clone_id ($db->get_all_Clone_id()) {
    print STDERR "\nDumping  clone  $clone_id\n";
    my $clone = $db->get_Clone($clone_id);
    foreach my $contig ($clone->get_all_Contigs()) {
	print STDERR "         contig ",$contig->id,"\n";
	print "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n\n";
	print "<!DOCTYPE game SYSTEM \"GAME.dtd\">\n\n";
	print "<game>\n";
	print " <annotation id=\"ENS:",$contig->id,"\n";
	print "   <name>ENS:",$contig->id,"</name>\n";
	print "    <version>1</version>\n";
	print "    <date> </date>\n";
	print " </annotation>\n";
	foreach my $gene ($contig->get_all_Genes()) {
	    print STDERR "         gene   ",$gene->id,"\n";
	    foreach my $trans ( $gene->each_Transcript ) {
		print STDERR"         trans. ",$trans->id,"\n";
		print " <feature_set annotation=\"",$trans->id,"\">\n";
		print "  <type>Transcript</type>\n";
		print "  <seq_relationship seq=\"",$contig->id,"\">\n";
		print "   <span>\n";
		my $switch = 1;
		my $trans_start;
		my $trans_end;
		foreach my $trans_exon ($trans->each_Exon) {
		    if ($trans_exon->contig_id eq $contig->id) {
			if ($switch == 1) {
			    $trans_start = $trans_exon->start;
			    $switch = 0;
			}
			$trans_end = $trans_exon->end;
		    }
		}
	        print "    <start>$trans_start</start>\n";
		print "    <end>$trans_end</end>\n";
		print "   </span>\n";
		print "  </seq_relationship>\n";
		print "  <feature_span>\n";
		print "  <type>translation</type>\n";
		print "   <seq_relationship seq=\"",$contig->id,"\">\n";
		print "    <span>\n";
		my $switch = 1;
		my $tr_start;
		my $tr_end;
		foreach my $tr_exon ($trans->translateable_exons) {
		    if ($tr_exon->contig_id eq $contig->id) {
			if ($switch == 1) {
			    $tr_start = $tr_exon->start;
			    $switch = 0;
			}
			$tr_end = $tr_exon->end;
		    }
		}
		print "     <start>$tr_start</start>\n";
		print "     <end>$tr_end</end>\n";
		print "    </span>\n";
		print "   </seq_relationship>\n";
		print "  </feature_span>\n";
		foreach my $exon ( $trans->each_Exon ) {
		    if ($exon->contig_id eq $contig->id) {
			print "  <feature_span>\n";
			print "  <type>Exon</type>\n";
			print "   <seq_relationship seq=\"",$exon->contig_id,"\">\n";
			print "    <span>\n";
			print "     <start>",$exon->start,"</start>\n";
			print "     <end>",$exon->end,"</end>\n";
			print "    </span>\n";
			print "   </seq_relationship>\n";
			print "   <evidence>\n";
			print "    <???>\n";
			print "   </evidence>\n";
			print "  </feature_span>\n";
		    }
		}
		print " </feature_set>\n";
		print "</game>\n";
	    }
	}
    }
}




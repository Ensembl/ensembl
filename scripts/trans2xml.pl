#!/usr/local/bin/perl

=head1 NAME

trans2xml

=head1 SYNOPSIS
 
 trans2xml.pl

=head1 DESCRIPTION

trans2xml produces XML file outputs of each contig, divded by transcripts

=cut
use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::SeqIO;

use Getopt::Long;

my $thost   = 'sol28';
my $tport   = '410000';
my $tdbname = 'ensdev';
my $format  = 'transcript';
my $user    = 'ensembl';

&GetOptions( 
	     'host:s'     => \$thost,
	     'port:n'     => \$tport,
	     'user:s'     => \$user,
	     'dbname:s'   => \$tdbname,
	     );

my $db = Bio::EnsEMBL::DBSQL::Obj->new( -user => $user, -db => $tdbname , -host => $thost );
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
		my $ex_cont_start;
		my $ex_cont_end;
		foreach my $tr_exon ($trans->translateable_exons) {
		    if ($tr_exon->contig_id eq $contig->id) {
			if ($switch == 1) {
			    $ex_cont_start = $tr_exon->start;
			    $switch = 0;
			}
			$ex_cont_end = $tr_exon->end;
		    }
		}
	        print "    <start>$ex_cont_start</start>\n";
		print "    <end>$ex_cont_end</end>\n";
		print "   </span>\n";
		print "  </seq_relationship>\n";
		print "  <feature_span>\n";
		print "  <type>translation</type>\n";
		print "   <seq_relationship seq=\"",$contig->id,"\">\n";
		print "    <span>\n";
		print "     <start>",$trans->translation->start,"</start>\n";
		print "     <end>",$trans->translation->end,"</end>\n";
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




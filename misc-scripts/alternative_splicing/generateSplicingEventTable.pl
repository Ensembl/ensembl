#!/usr/local/ensembl/bin/perl -w

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

 generateSplicingEventTable.pl

=head1 AUTHORS

 Gautier Koscielny (koscieln@ebi.ac.uk)

=head1 SYNOPSIS

  generateSplicingEventTable.pl -dbhost host -dbuser ensro -dbname homo_sapiens_core_62_37g

=head1 DESCRIPTION

This script generates a tabular separated table containing all the splicing events for a set of genes.

here is an example commandline

  ./generateSplicingEventTable.pl -dbhost host -dbuser user -dbname my_db -dbpass ****

=head1 OPTIONS

    -dbhost         host name for database (gets put as host= in locator)
		-dbport         port (gets put as port= in locator)
    -dbname         what name to connect to (dbname= in locator)
    -dbuser         what username to connect as (dbuser= in locator)
    -dbpass         what password to use (dbpass= in locator)
    -help           displays this documentation with PERLDOC

=cut

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SplicingEventAdaptor;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::AttributeAdaptor;
use Bio::EnsEMBL::DBSQL::ExonAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);


my $host   = '';
my $port   = '3306';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my $file   = undef;
my $help;
my @coord_system;

&GetOptions(
		'dbhost:s'   => \$host,
		'dbport:n'   => \$port,
		'dbname:s'   => \$dbname,
		'dbuser:s'   => \$dbuser,
		'file:s'      => \$file,
		'dbpass:s'   => \$dbpass,
		'h|help'     => \$help,
		) or ($help = 1);

if(!$host || !$dbuser || !$dbname) { 
		print STDERR "Can't get any information without database details\n";
		print STDERR "-dbhost $host -dbuser $dbuser -dbname $dbname ".
				" -dbpass $dbpass\n";
		$help = 1;
}

if ($help) {
		exec('perldoc', $0);
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
		(-dbname => $dbname,
		 -host   => $host,
		 -user   => $dbuser,
		 -port   => $port,
		 -pass   => $dbpass);

my $gene_adaptor = $db->get_GeneAdaptor();
my $se_adaptor = $db->get_SplicingEventAdaptor();
my $attrib_adaptor = $db->get_AttributeAdaptor();
my $exon_adaptor = $db->get_ExonAdaptor();
my @stable_gene_ids = @{$gene_adaptor->list_stable_ids()};
my $size = scalar @stable_gene_ids;
print STDERR "Number of stable ids:\t" . $size . "\n";

my $count = 0;

# example genes 
my @genes = ("ENSG00000000419", "ENSG00000006282");

# don't know how to get this information easyly from the core API

my $attrib_names = {
		II => "Intron isoform",
		IR => "Intron retention",
		A5SS => "Alternative 5' splice site",
		A3SS => "Alternative 3' splice site",
		CE => "Cassette exon",
		MXE => "Mutually exclusive exons",
		EI => "Exon isoform",
		ALE => "Alternative last exon",
		AFE => "Alternative first exon"
				
};

print " AS name\tType\tShort desc\tlocation\tdetails\n";
print "+--------------------------------------------------------------------------+\n";
for my $id (@genes) {
		
		$count++;
		
		my $gene = $gene_adaptor->fetch_by_stable_id($id);
		
		my $gene_id = $gene->display_id();
		my $biotype = $gene->biotype();
		my $chr = $gene->slice->seq_region_name();
		my $strand = $gene->strand();
		my $start = $gene->start();
		my $end = $gene->end();
		
		# get all splicing events
		
		my @irs = ();

		for my $se ( sort { if ($a->type() eq $b->type()) { if ($a->strand == 1) { $a->start() <=> $b->start() } else { $a->end() <=> $b->end() } } } @ { $se_adaptor->fetch_all_by_Gene($gene) }) {
				
				next if ($se->type() eq 'CNE');

				my $location = ($se->strand() == 1) ? ("+\t$chr:" . $se->start() . "-" . $se->end()) : ("-\t$chr:" . $se->end() . "-" . $se->start());
				
				if ($se->type() eq 'IR') {
						# special case for IR
						push @irs, $se;
				} else {
						# display of the exons depends on the event
						print $se->type() . "\t" . $attrib_names->{$se->type()} . "\t$location";
						print "\t" . get_event_display($se) . "\n";
				}
				
		}

		# finally, IR display if any
		print get_IR_display(\@irs, $chr) if (scalar(@irs) > 0);
		
		# last if ($count >=3);
}


sub get_IR_display {
		my ($ref, $chr) = @_;

		my $s = '';
		my $exons = {};

		for my $se (sort { if ($a->strand == 1) { $a->start() <=> $b->start() } else { $a->end() <=> $b->end() } } @{ $ref }) {
				
				my @features = @{ $se->get_all_Features() };
				my %f_map = map { $_->exon_id() => $_ } @features;
				
				foreach my $exon_id (sort { $f_map{$a}->feature_order() <=> $f_map{$b}->feature_order() } keys %f_map) {
						
						my $f = $f_map{$exon_id};
						my $exon = $exon_adaptor->fetch_by_dbID($exon_id);
						my $retained = $exon->stable_id() . "(" . $f->start() . "-" . $f->end() . ")";

						if ($f->feature_order() == 1) {
								
								if (!defined($exons->{$exon->stable_id()})) {
										my $location = ($se->strand() == 1) ? ("+\t$chr:" . $f->start() . "-" . $f->end()) : ("-\t$chr:" . $f->end() . "-" . $f->start());
										$s .= $se->type() . "\t" . $attrib_names->{$se->type()} . "\t$location\tIntron retained in: $retained\n";
										$exons->{$exon->stable_id()} = {
												text => $se->type() . "\t" . $attrib_names->{$se->type()} . "\t$location\tIntron retained in: $retained\n",
												start => $f->start(),
												end => $f->end()
										}
								}
								
						}
				}
				
		}
		
		return $s;
}

sub get_event_display {
		my $se = shift;
		
		my $s = '';
		my @features = @{ $se->get_all_Features() };
		my %f_map = map { $_->exon_id() => $_ } @features;
		#print join("...", keys %f_map) . "\n";
		
		if ($se->type() eq 'CE') {
				
				# what are the exons skipped and their coordinates?
				$s = (scalar(keys %f_map) > 1) ? "Skipped exons: " : "Skipped exon: ";
				my $nb = 0;
				my @sl = ();
				foreach my $exon_id (sort { $f_map{$a}->feature_order() <=> $f_map{$b}->feature_order() } keys %f_map) {
						my $f = $f_map{$exon_id};
						my $exon = $exon_adaptor->fetch_by_dbID($exon_id);
						push @sl, $exon->stable_id() . "(" . $f->start() . "-" . $f->end() . ")";
						
				}
				$s .= join("; ", @sl);
				
		} elsif ($se->type() eq 'MXE' || $se->type() eq 'AFE' || $se->type() eq 'ALE' ) {

				my $nb = 0;
				my @s1 = ();
				my @s2 = ();

				foreach my $exon_id (sort { $f_map{$a}->feature_order() <=> $f_map{$b}->feature_order() } keys %f_map) {
						my $f = $f_map{$exon_id};
						my $exon = $exon_adaptor->fetch_by_dbID($exon_id);
						my $s = $exon->stable_id() . "(" . $f->start() . "-" . $f->end() . ")";
						my $v = ($se->type() eq 'MXE') ? $f->transcript_association() : $f->feature_order();
						if ($v == 1) {
								push @s1, $s;
						} else {
								push @s2, $s;
						}
						
				}
				$s = "Alternative exons: " . join("; ", @s1) . " / " . join("; ", @s2);
				
		} elsif  ($se->type() eq 'II') {

				my $nb = 0;
				my @s1 = ();
				my @s2 = ();

				foreach my $exon_id (sort { $f_map{$a}->feature_order() <=> $f_map{$b}->feature_order() } keys %f_map) {
						my $f = $f_map{$exon_id};
						my $exon = $exon_adaptor->fetch_by_dbID($exon_id);
						my $s = $exon->stable_id() . "(" . $f->start() . "-" . $f->end() . ")";
						if ($f->transcript_association() == 1) {
								push @s1, $s;
						} else {
								push @s2, $s;
						}
						
				}

				$s = "Flanking exons: " . join("; ", @s1) . " / " . join("; ", @s2); 

		} elsif  ($se->type() eq 'EI') {
				
				my $nb = 0;
				my @s1 = ();
				my @s2 = ();

				foreach my $exon_id (sort { $f_map{$a}->feature_order() <=> $f_map{$b}->feature_order() } keys %f_map) {
						my $f = $f_map{$exon_id};
						my $exon = $exon_adaptor->fetch_by_dbID($exon_id);
						my $s = $exon->stable_id() . "(" . $f->start() . "-" . $f->end() . ")";
						if ($f->feature_order() == 1) {
								push @s1, $s;
						} else {
								push @s2, $s;
						}
						
				}

				$s = "Alternative exons: " . join("; ", @s1) . " / " . join("; ", @s2);

		} elsif  ($se->type() eq 'A5SS' ||  $se->type() eq 'A3SS') {

				my $nb = 0;
				my @s1 = ();
				my @s2 = ();
				
				foreach my $exon_id (sort { $f_map{$a}->feature_order() <=> $f_map{$b}->feature_order() } keys %f_map) {
						
						my $f = $f_map{$exon_id};
						my $exon = $exon_adaptor->fetch_by_dbID($exon_id);
						my $s = $exon->stable_id() . "(" . $f->start() . "-" . $f->end() . ")";

						if ($f->feature_order() == 1) {
								push @s1, $s;
						} else {
								push @s2, $s;
						}
						
				}

				$s = "Alternative flanking exons: " . join("; ", @s1) . " / " . join("; ", @s2);

		}

		return $s;

}


exit 0;

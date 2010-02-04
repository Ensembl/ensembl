#!/usr/local/ensembl/bin/perl -w

=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

 Fetch_gff.pl

=head1 AUTHORS

 Gautier Koscielny (koscieln@ebi.ac.uk)

=head1 SYNOPSIS

  Fetch_gff.pl -dbhost host -dbuser ensro -dbname homo_sapiens_core_57_37b

=head1 DESCRIPTION

This script generates a GFF file containing all the transcript variants 
from an Ensembl core database.

here is an example commandline

  ./Fetch_gff.pl -dbhost host -dbuser user -dbname my_db -dbpass ****

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)
    -dbname    what name to connect to (dbname= in locator)
    -dbuser    what username to connect as (dbuser= in locator)
    -dbpass    what password to use (dbpass= in locator)
    -chr       which chromosome (optional)
    -help      displays this documentation with PERLDOC

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use POSIX qw(ceil);

use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

{ # block to avoid namespace pollution
  my $host   = '';
  my $port   = '';
  my $dbname = '';
  my $dbuser = '';
  my $dbpass = '';
	my $chr    = undef;
  my $help;
  my @coord_system;

  &GetOptions(
              'dbhost:s'   => \$host,
              'dbport:n'   => \$port,
              'dbname:s'   => \$dbname,
              'dbuser:s'   => \$dbuser,
							'chr:s'      => \$chr,
              'dbpass:s'   => \$dbpass,
              'h|help'     => \$help,
             ) or ($help = 1);

  if(!$host || !$dbuser || !$dbname || !$dbpass){
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
	my @stable_gene_ids = undef;
	my $size = 0;

	if (defined($chr)) {

			my $slice_adaptor = $db->get_SliceAdaptor();
			my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
			@stable_gene_ids = @{ $gene_adaptor->fetch_all_by_Slice($slice) };
			$size = scalar @stable_gene_ids;
			print STDERR "Number of stable genes on region $chr:\t" . $size . "\n";

	} else {

			@stable_gene_ids = @{$gene_adaptor->list_stable_ids()};
			$size = scalar @stable_gene_ids;
			print STDERR "Number of stable ids:\t" . $size . "\n";
	}

	for my $id (@stable_gene_ids) {
			
			my $gene = ($chr) ? $id : $gene_adaptor->fetch_by_stable_id($id);

			my $gene_id = $gene->display_id();
			my $biotype = $gene->biotype();
			my $chr = $gene->slice->seq_region_name();
			my $strand = $gene->strand();
			my $start = $gene->start();
			my $end = $gene->end();
			
			my @transcripts = @{$gene->get_all_Transcripts()};
			for my $transcript (@transcripts) {
					
					my $transcr_id = $transcript->display_id() ; ;
					
					#Get the exons + print info.
					my $exons = $transcript->get_all_Exons() ;
					
					foreach my $exon (@$exons) {
							my $exon_id = $exon->display_id() ;
							my $exon_start = $exon->start() ;
							my $exon_end = $exon->end() ;
							my $exon_std = $exon->strand() ;
							my $slice = $exon->slice->seq_region_name();
							$exon_std =~ s/-1/-/ ;
							$exon_std =~ s/1/+/ ;
							print "$chr\tEnsembl\texon\t$exon_start\t$exon_end\t.\t$exon_std\t.\tgene_id \"$gene_id\"; transcript_id \"$transcr_id\"; exon_id \"$exon_id\"\n" ;
					}
			}
			
	}
	
	exit 0;

}

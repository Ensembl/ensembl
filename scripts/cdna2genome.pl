#!/usr/local/bin/perl

use strict;

#  perl ../scripts/cdna2genome.pl -timdb ~/th/unfinished_ana/sanger/data/dJ718J7/dJ718J7.seq

=head1 NAME

cdna2genome

=head1 SYNOPSIS

    cdna2genome fasta-file-genome list-of-accession-numbers

    cdna2genome -timdb fasta-file-genome-sequence (*)

    -output [gff/embl]
    -nosim4 - do not show sim4 exons/exonsets
    -nocds  - do not show cds 
    -timdb  - switch to turn on ensembl pipeline processing

=head1 DESCRIPTION

This module compares a genomic sequence as fasta file
to a series of cDNAs designated by EMBL accession numbers.
These EMBL accession numbers are extracted and then
aligned using Sim4. Then the EMBL based CDS line is used
to provide coding sequence prediction.

The results are then merged to remove completely identical 
predictions and given out in GFF format

(*) full file specification of a sequence file containing multiple
contigs for an ensembl pipeline analysis directory.  Each contig is
processed in turn and the accessions of vertrna hits to these clones
are pulled out by scanning the hits to each predicted peptide for that
contig.  Output written to file in analysis directory
contig.sim4_vert.gff

=cut

use Bio::Tools::Sim4::Results;
use Bio::AnnSeq;
use Bio::SeqIO;
use Bio::AnnSeqIO;

use Bio::EnsEMBL::Analysis::Log;

use Getopt::Long;


my $log = new Bio::EnsEMBL::Analysis::Log;
$log->dump_to_error(1);


my $output = 'gff';
my $nosim4 = 0;
my $nocds  = 0;
my $nologerr = 0;
my $timdb;

&GetOptions( 
	     'output:s' => \$output,
	     'nosim4'   => \$nosim4, 
	     'nocds'    => \$nocds,
             'nolog'    => \$nologerr,
	     'timdb'    => \$timdb,
	     );


my $format = 'Fasta';
my $file  = shift;

my %top_sf;

$file || do { system("perldoc $0"); exit(0) };

# get directory and clonename if timdb
my($dir,$clone,@files);
if($timdb){
    if($file=~/(.*)\/([^\/]+)\.seq/){
	$dir=$1;
	$clone=$2;
	local *TMPDIR;
	opendir(TMPDIR,$dir);
	@files=readdir(TMPDIR);
	closedir(TMPDIR);
    }else{
	die "cannot parse $file as timdb seq file";
    }
}

my $seqin = Bio::SeqIO->new( -format => $format , -file => $file );

# file could have more than sequences, so process them in turn
my $seq;
my $nseq;
while($seq = $seqin->next_seq()){
    $nseq++;

    my @embl_acc;
    my $OUTFILE;
    if($timdb){
	# get list of targets from an unfinished_ana directory if specified

	# get contig name and output file
	my $contig=$seq->id;
	#print "$contig\n";

	# loop over peptide hit files for this contig
	my %acc;
	foreach my $file (@files){
	    if($file=~/$contig\.(\d+)\.tblastn_vert\.msptmp$/){
		#print "  $file\n";
		my $pepid=$1;
		local *TMPFILE;
		open(TMPFILE,"$dir/$file") || die "cannot open $file";
		while(<TMPFILE>){
		    s/^\s+//;
		    my @fields=split(/\s+/,$_);
		    $acc{$fields[8]}.=$pepid.",";
		}
		close(TMPFILE);
	    }
	}
	@embl_acc=(keys %acc);
	
	# next unless msptmp hits
	next unless scalar(@embl_acc);

	#my $outfile="$dir/$contig.sim4_vert.$output";
	#open(SIM4OUT,">$outfile") || die "cannot open $outfile";
	#$OUTFILE=\*SIM4OUT;

	$OUTFILE = \*STDOUT;

	# debug
	#foreach my $acc (@embl_acc){
	#    print $OUTFILE "$contig $acc $acc{$acc}\n";
	#}
 
    }else{
	@embl_acc=@ARGV;
	$OUTFILE=\*STDOUT;
    }

    &_do_cdna2genome($OUTFILE,\@embl_acc);

    if($timdb){
	close($OUTFILE);
    }

}

# if never found a sequence, report
unless($nseq){
    die "no sequence read from $file";
}

# this is the end of the script now
exit 0;

sub _do_cdna2genome{
    my($OUTFILE,$rembl_acc)=@_;

    my $aseq = Bio::AnnSeq->new();
    $aseq->seq($seq);

    #my $OUTFILE;

  ACC: 
    foreach my $embl_acc ( @$rembl_acc ) {

	# get out the cdna.
	my $aseqio = Bio::AnnSeqIO->new( -format => 'EMBL', 
					 -file => "efetch -a em:$embl_acc |");    
	my $cdna;
	eval {
	    $cdna = $aseqio->next_annseq();
	};

	# skip if we can't load it or no cdna.
	if( $@ || !$cdna) {
	    $log->store("Unable to read $embl_acc. Skipping. (Apologies)\n$@\n");
	    next ACC;
	}

	# get out CDS coordinates.
	
	my $cds_start;
	my $cds_end;
	
	foreach my $sf ( $cdna->top_SeqFeatures ) {
	    if( $sf->primary_tag() eq 'CDS' ) {
		$cds_start = $sf->start;
		$cds_end   = $sf->end;
	    }
	}
	
	if( !defined $cds_start ) {
	    print STDERR "EMBL $embl_acc has not single CDS line. Skipping\n";
	}
	
    
	# wrap everything in an eval to trap exceptions
	eval {
	    # assumme file is in fasta format for sim4 stuff.
	    
	    my $seqout = Bio::SeqIO->new( -format => 'Fasta', -file => ">/tmp/sim4.cdna.$$");
	    $seqout->write_seq($cdna->seq);
	    $seqout = undef;
	    print STDERR "file $file vs /tmp/sim4.cdna.$$ \n";

	    my $sim4 = Bio::Tools::Sim4::Results->new( -file => "sim4 $file /tmp/sim4.cdna.$$ |" );
	    

	    foreach my $es ( $sim4->each_ExonSet ) {
		
		# add this exon onto the AnnSeq anyway.
		
		$aseq->add_SeqFeature($es);
		
		# we need to process this into a set of CDS regions. 
		# we need to get the first exon.
		
		my $seen_start = 0;
		my @exons = $es->each_Exon;

		if( $#exons == -1 ) {
		    die("No exons in $embl_acc sim4 run!");
		}

		
		# sort by whether homology is start/end.
		
		@exons = sort { $a->homol_SeqFeature->start <=> $b->homol_SeqFeature->end } @exons;
		
		my $strand;
		
		$strand = $exons[0]->homol_SeqFeature->strand();
		
		my $exon; # temporay variable of the current exon
		while( $exon = shift @exons ) {
		    if( $exon->homol_SeqFeature->start < $cds_start && 
			$exon->homol_SeqFeature->end > $cds_start ) {
			$seen_start =1;
			last;
		    }
		}

		if( $seen_start == 0 ) {
		    $log->store("Have not seen start. Skipping an exon set");
		    next;
		}
		
	    
		#
		# Build a top level Generic SeqFeature object
		# and add in the individual exons
		#
	    
		my $top_sf = Bio::SeqFeature::Generic->new();
		$top_sf->primary_tag('CDS_span');
		$top_sf->source_tag('Sim4_CDS');
		$top_sf->_parse->{'parent_homogenous'} =1;
		$top_sf->strand($strand);
		$top_sf->add_tag_value('match_cdna',$embl_acc);

		# first exon. Has to change start.
		my $sf = Bio::SeqFeature::Generic->new();
		if( $strand == 1 ) {
		    $sf->start($exon->start + ($cds_start - $exon->homol_SeqFeature->start));
		    $sf->end($exon->end);
		} else {
		    $sf->start($exon->start);
		    $sf->end($exon->end - ($cds_start - $exon->homol_SeqFeature->start));
		}
		
		$sf->strand($strand);
		$sf->primary_tag('CDS');
		$sf->source_tag('Sim4_CDS');
		$sf->add_tag_value('match_cdna',$embl_acc);
		$top_sf->add_sub_SeqFeature($sf,'EXPAND');
		
		my $seen_end = 0;
		my $end_exon;
		foreach $exon ( @exons ) {
		    # the remainder
		    if( $exon->homol_SeqFeature->start() < $cds_end && 
			$exon->homol_SeqFeature->end() > $cds_end ) {
			$seen_end = 1;
			$end_exon = $exon;
			last;
		    }

		    my $sf = Bio::SeqFeature::Generic->new();
		    if( $strand == 1 ) {
			$sf->start($exon->start);
			$sf->end($exon->end);
		    } else {
			$sf->start($exon->start);
			$sf->end($exon->end);
		    }

		    $sf->strand($strand);
		    $sf->primary_tag('CDS');
		    $sf->source_tag('Sim4_CDS');
		    $sf->add_tag_value('match_cdna',$embl_acc);
		    $top_sf->add_sub_SeqFeature($sf,'EXPAND');
		}

		$exon= $end_exon;

		# last exon.
		if( $seen_end == 0 ) {
		    $log->store("Did not see end exon. Yuk.\n");
		} else {
		    
		    my $sf = Bio::SeqFeature::Generic->new();
		    if( $strand == 1 ) {
			$sf->start($exon->start);
			$sf->end($exon->start + ($cds_end - $exon->homol_SeqFeature->start));
		    } else {
			$sf->start($exon->end - ($cds_end - $exon->homol_SeqFeature->start));
			$sf->end($exon->end);
		    }
		    
		    $sf->strand($strand);
		    $sf->primary_tag('CDS');
		    $sf->source_tag('Sim4_CDS');
		    $sf->add_tag_value('match_cdna',$embl_acc);
		    $top_sf->add_sub_SeqFeature($sf,'EXPAND');
		}

		# add in this feature
		
		$top_sf->strand($strand);
		$top_sf{$top_sf} = $top_sf;
		
	    }
	};

	# catch exceptions
	
	if( $@ ) {
	    $log->store("unable to process $embl_acc\n$@\n");
	}
	
    }

    #
    # Sort by size
    #

    my @size = sort { $a->length <=> $b->length } values %top_sf;

    #
    # Go down the list, marking off guys which are subsummed by deleting from the hash
    #
    
    foreach my $sf ( @size ) {
	foreach my $test ( values %top_sf ) {
	    if( $test == $sf ) {
		next;
	    }
	    
	    if( compare_top_seqfeature($test,$sf) == 1 ) {
		# remove sf
		delete $top_sf{$sf};
	    }
	}
    }

    #
    # Add them to AnnSeq object
    #

    foreach my $sf ( values %top_sf ) {
	$aseq->add_SeqFeature($sf);
    }
    

    if( $output eq 'embl' ) {
	my $aout = Bio::AnnSeqIO->new( -format => 'EMBL', 
				       -fh => $OUTFILE );
	$aout->write_annseq($aseq);
    } else {
	foreach my $sf ( $aseq->all_SeqFeatures() ) {
	    if( $nosim4 ) {
		if( $sf->source_tag() eq 'Sim4' ) { next; }
	    }
	    if( $nocds ) {
		if( $sf->source_tag() eq 'Sim4_CDS' ) { next; }
	    }
	    
	    print $OUTFILE $sf->gff_string(), "\n";
	}
    }
    
    if( $nologerr == 0 ) {
	print STDERR "Messages logged during executation. Switch off with -nolog\n";
	foreach my $mess ( $log->each_LogMessage() ) {
	    print "Message: ",$mess->text," by: ",$mess->author, " at:\n", $mess->stack, "\n\n";
	}
    }

}

##################
# SubRoutines
##################

sub compare_top_seqfeature {
    my ($one,$two) = @_;

    # is two completely contained by one

    if( $one->strand ne $two->strand ) {
	return -1;
    }
    my $strand = $one->strand;

    if( $one->start > $two->start || $one->end < $two->end ) {
	# can't be contained
	return -1;
    }
    
    # first and last sub features handled differently.

    my @exons = $two->sub_SeqFeature;
    my $first = shift @exons;
    my $last  = pop @exons;

    my @matched_exons = $one->sub_SeqFeature;

    #
    # Order is important. It is not good enough just that
    # exons in two match one, as if one has an additional exon,
    # then they are separate
    #

    while ( my $current_exon = shift @matched_exons ) {
	if( $strand == 1 ) {
	    if( $current_exon->start < $first->start && $current_exon->end > $first->start ) {
		if( $first->end != $current_exon->end ) {
		    return -1;
		} else {
		    last;
		}
	    }
	} else {
	    if( $current_exon->end > $first->end && $current_exon->start < $first->end ) {
		if( $first->start != $current_exon->start ) {
		    return -1;
		} else {
		    last;
		}
	    }
	}
    }

    if( $#matched_exons == -1 ) {
	return -1;
    }

    #
    # ok. Now each exon should agree 
    #

    while( my $current_two_exon = shift @exons ) {
	my $current_exon = shift @matched_exons;
	if( !defined $current_exon ) {
	    return -1;
	}

	if( $current_two_exon->start != $current_exon->start ||
	    $current_two_exon->end != $current_exon->end ) {
	    return -1;
	}
    }

    #
    # final exon
    #

    my $last_one_exon = shift @matched_exons;

    if( $strand == 1 ) {
	if( $last->start != $last_one_exon->start || $last->end > $last_one_exon->end )  {
	    return -1;
	}
    } else {
	if( $last->end != $last_one_exon->end || $last->start < $last_one_exon->start )  {
	    return -1;
	}
    }

    return 1;
}
	    


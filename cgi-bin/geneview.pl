#!/usr/local/bin/perl 

# makes GFF stuff for a contig.

BEGIN {
    push(@INC,"../modules");
    push(@INC,"../../../bioperl-live");
}

use CGI;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::Gene;
use strict;

my $q = new CGI;
print $q->header();

my $geneid = $q->param('gene');
my $gene;
my $db;

printHeader($geneid);

eval {
    $db = new Bio::EnsEMBL::DBSQL::Obj( -user   => 'ensro', 
				       -dbname => 'ensembl' , 
				       -host   => 'localhost');


    $gene = $db->get_Gene($geneid);
};

if( $@ ) {
    print(STDERR "Error $@\n");
} else {



    foreach my $transcript ($gene->each_Transcript) {

	$db->get_supporting_evidence($transcript->each_Exon);

	print("<h2>Transcript " . $transcript->id . "</h2>\n");

	print("Number of exons in transcript = " . scalar($transcript->each_Exon) . "<br>\n");

	print("<h3>Exons are :</h3>\n");

	my $count = 1;
	print("<ul>\n");
	foreach my $exon ($transcript->each_Exon) {
	    print("<li>Exon $count : " . $exon->id . "</li>\n");
	    $count++;
	}
	print("</ul>\n");
	print("<h3>Supporting evidence for exons</h3>\n");


	print("Below is a table showing which database hits have overlaps\n");
	print("with each exon in the transcript.  The database hits are the results of a series of blast runs against genscan predicted peptides.<p>\n");

	    

	# Hash the evidence by feature hid
	my %hid;
	my %type;
	my @features;

	foreach my $exon ($transcript->each_Exon) {
	    foreach my $f ($exon->each_Supporting_Feature) {
		if (!defined($hid{$f->hseqname})) {
		    $hid{$f->hseqname} = [];
		}
		$type{$f->hseqname} = $f->analysis->db;
		push(@{$hid{$f->hseqname}},$exon);
		push(@features,$f);
	    }
	}

	# Now do some sorting on the features;
	@features = sort { $a->hseqname <=> $b->hseqname ||
			       $a->analysis->db <= $b->analysis->db ||
			       $#{$hid{$b->hseqname}} <=> $#{$hid{$a->hseqname}} } @features;

	# now rank the hid
	my $prev;
	my @hid;

	foreach my $f (@features) {
	    if (!$prev || ($f->hseqname ne $prev)) {
		push(@hid,$f->hseqname);
		$prev = $f->hseqname;
	    }
	}

	print("<pre>\n");
	# First the top line for the exons
	printf("%25s","Exon number =>");
	printf("%15s","");

	my $count = 1;

	foreach my $exon ($transcript->each_Exon) {
	    printf("%8s",$count);
	    $count++;
	}

	print("\n");

	foreach my $hid (@hid) {
	    printf("%-25s",$hid);
	    printf("%-15s",$type{$hid});

	    foreach my $exon ($transcript->each_Exon) {
		my $feature;

		foreach my $f ($exon->each_Supporting_Feature) {
		    if ($f->hseqname eq $hid) {
			$feature = $f;
		    }
		}

		if (defined($feature)) {
		    printf("%8s","=======");
		} else {
		    printf("%8s"," ");
		}
	    }
	    print("\n");
	}
	print("</pre>\n");
    }

}

sub printHeader {
    my ($geneid) = @_;
    print("<h1>Supporting Evidence for gene $geneid</h1>\n");
    print("The EnsEMBL database is now being queried to find the supporting evidence for the gene $geneid.  Please wait as this could take several seconds<p>\n");
}

    

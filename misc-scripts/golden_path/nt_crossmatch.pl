#!/usr/local/ensembl/bin/perl

#This script reads a file which has nt contig and clone coordinates
#then reads a fasta file of all nt contigs
#then runs crossmatches between NT contigs and each raw contig mapped to be within that NT

use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;

$|=1;

#This is the NT-clone agp file
my $ctg_coords=shift(@ARGV);

#This is the NT fasta file
my $ntfile=shift(@ARGV);
$ntfile =~ /(\S+)\.fa/;
my $nt= $1;

my $gp="/work2/elia/oct07/data/$nt.agp";


my $seqin  = Bio::SeqIO->new(-file => "/work2/elia/oct07/data/$ntfile" , '-format' => 'Fasta');
print STDERR "Reading NT fasta file...\n";

my $ntseq = $seqin->next_seq();

my %nt_clone;

open (CTG,"<$ctg_coords");
print STDERR "Reading ctg coordinate file...\n";
while (<CTG>) {
    my @tabs=split(/\t/);
    push (@{$nt_clone{$tabs[0]}},$tabs[6]);
}
close (CTG);

#Hash holding the seq object for each NT contig
my %nt_hsp;
my %nt_fasta;

$ntseq->id =~ /ref\|(\S+)\|(\S+)/;
my $id=$1;
my $hsp=$2;
$nt_hsp{$id}=$hsp;
$ntseq->display_id($id);

my $db=Bio::EnsEMBL::DBSQL::Obj->new(-dbname=>'arne_oct07_tim',-host=>'ecs1b',-user=>'lijnzaad');
print STDERR "Starting crossmatches...\n";

open (GP,">$gp");

print STDERR "Going to crossmatch $nt:\n";
my @matches;
foreach my $cl (@{$nt_clone{$nt}}) {
    print STDERR "Getting $cl from db\n";
    $cl =~ /(\S+)\.(\d+)/;
    my $clone=$db->get_Clone($1);
    if (($clone->embl_version != $2) && ($clone->embl_version !=0)) {
	print STDERR "BUG: $cl has version $2 in the flatfile and ".$clone->embl_version." in the database\n";
	next;
    }
    print STDERR "Going to crossmatch $cl with $nt...\n";
    my $c=0;
    my @fp;
    my %offset;
    print STDERR "Getting here...\n";
    foreach my $contig ($clone->get_all_Contigs) {
	$c++;
	$offset{$contig->id}=$contig->embl_offset;
	my $seq=$contig->primary_seq->seq();
	$seq = uc ($seq);
	$seq =~ s/[^ATCG]/N/g;
	my $raw_seq = Bio::PrimarySeq->new( -display_id => $contig->primary_seq->id, -seq => $seq);
	my $nts=$ntseq->seq;
	$nts =~ s/[^ATCG]/N/g;
	$nts = uc ($nts);
	$ntseq->seq($nts);
	print STDERR "Crossmatching ".$contig->id." with $nt\n";
	my $score=200;
	my $cross = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new( 
								       -nocopy => 1,
								       -seq1 => $raw_seq,
								       -seq2 => $ntseq,
								       -score => $score
								       );
	print STDERR "Using score: $score\n";

	$cross->run();
        print STDERR "Finished crossmatch\n";
	push (@fp,$cross->output);
    }
    if ($c > 1) {
	print STDERR "BUG: $cl has $c contigs even though it is finished\n";

    }
    my $size=scalar(@fp);
    if ($size == 0) {
	print STDERR "BUG: $cl-$nt crossmatch did not return any feature pairs, need to rerun with lowers threshold";
    }
    elsif ($size == 1) {
	#Easy case, one feature pair between the two
	print STDERR "One feature pair for $nt-$cl\n";
	my $f = $fp[0];

	my $cv= $clone->id.".".$clone->embl_version;
	my $match;
	$match->{clone} = $cv;
	$match->{NT} = $nt;
	$match->{ntstart} = $f->hstart;
	$match->{ntend} = $f->hend;
	my $cl_start = $f->start+$offset{$f->seqname}-1; #Adjust for offset!
	$match->{clstart} = $cl_start;
	my $cl_end = $f->end+$offset{$f->seqname}-1; #Adjust for offset!
	$match->{clend} = $cl_end;
	if ($f->hstrand == 1) {
	    #print STDERR "...forward...\n";
	    $match->{strand} = '+1';
	}
	else {
	    #print STDERR "...reverse...\n";
	    $match->{strand} = '-1';
	}
	push (@matches,$match);
    }
    else { 
	print STDERR "More than one feature pair ($size) for $cl-$nt, getting longest\n";

	my $cv= $clone->id.".".$clone->embl_version;

	#Take longest crossmatch!
	my $length=0;

	#This loop finds longest match
	foreach my $f (@fp) {
	    my $match;
	    $match->{clone} = $cv;
	    $match->{NT} = $nt;
	    print STDERR "Checking feature on nt ".$f->hstart." ".$f->hend." and CLONE ".$f->start." ".$f->end."\n";
	    #if (($f->hend-$f->hstart)> $length) {
	    $match->{ntstart} = $f->hstart;
	    $match->{ntend} = $f->hend;
	    my $cl_start = $f->start+$offset{$f->seqname}-1;
	    $match->{clstart} = $cl_start;
	    my $cl_end = $f->end+$offset{$f->seqname}-1;
	    $match->{clend} = $cl_end;
	    if ($f->hstrand == 1) {
		$match->{strand} = '+1';
	    }
	    else {
		$match->{strand} = '-1';
		
	    }
	    #$length=$f->hend-$f->hstart;
	    #}
	    push (@matches,$match);
	}
    }
}
my $size = scalar (@matches);
print STDERR "Got $size matches for this NT contig\n";

my @validated;
my %internal;
if ($size>1) {

    #Start by sorting matches by start in nt
    my @sorted = sort {$a->{ntstart} <=> $b->{ntstart}} @matches;

    #This loop takes care of removin internal matches (i.e. matches within matches

    #We shift out one by one the matches from the array, and compare
    #them to all other subsequenet matches (ordered by ntstart)
    while (my $match = shift @sorted) {
	my $key1=$match->{clone}."-".$match->{clstart};
	if (exists $internal{$match->{$key1}}) {
	    #print STDERR "Skipping match with clone ".$match->{clone}.", because already detected as internal ".$match->{clstart}."-".$match->{clend}."\n";
	    next;
	}
	#print STDERR "Checking ".$match->{NT}."\t".$match->{ntstart}."\t".$match->{ntend}."\t".$match->{clone}."\t".$match->{clstart}."\t".$match->{clend}."\t".$match->{strand}."\n";

	#Check if subsequent  matches are contained in this match 
	foreach my $other (@sorted) {
	    my $key2=$other->{clone}."-".$other->{clstart};
	    #Check if contaied and not set to internal already
	    if (($other->{ntend} < $match->{ntend}) && (! exists $internal{$key2})) {
		#Store internal matches in %internal
		$internal{$key2}=1;
		#print STDERR "Skipping ".$other->{NT}."\t".$other->{ntstart}."\t".$other->{ntend}."\t".$other->{clone}."\t".$other->{clstart}."\t".$other->{clend}."\t".$other->{strand}."\n";
	    }
	}
	
	#All matches that are not internal, go into validated arrau
	if (! exists $internal{$key1}) {
	    push @validated,$match;
	}
    }
}
else { 
    push (@validated,@matches);
}

#Now we have a set of non-redundant matches sorted by nt_start
#Now we have to go through them to assign golden path switch points
#We write the determined golden path to a file, and try to detect gaps if any

print STDERR "\nFinal set:\n";
foreach my $match (@validated) {
    print STDERR $match->{NT}."\t".$match->{ntstart}."\t".$match->{ntend}."\t".$match->{clone}."\t".$match->{clstart}."\t".$match->{clend}."\t".$match->{strand}."\n";
   

}

my $size=scalar(@validated);
my $hsp=$nt_hsp{$nt};
print STDERR "The golden path will have $size elements for NT contig $nt:\n";
my ($ntstart,$ntend,$cl,$clstart,$clend,$strand);

if ($size == 1) {
    #Write golden path directly from match results
    my $match=$validated[0];
    $ntstart = $match->{ntstart};
    $ntend = $match->{ntend};
    $cl = $match->{clone};
    $clstart = $match->{clstart};
    $clend = $match->{clend};
    $strand = $match->{strand};
    print STDERR "$nt\t$hsp\t$ntstart\t$ntend\t$strand\tF\t$cl\t$clstart\t$clend\n";
    print GP "$nt\t$hsp\t$ntstart\t$ntend\t$strand\tF\t$cl\t$clstart\t$clend\n";
}
else {
    #Loop through this final validated set to deduce ("walk trough") golden path
    for (my $i=0; $i<$size; $i++) {
	#print STDERR "looking at $i size...\n";

	#Variables to be written to golden path file
	my $match=$validated[$i];
	my $next=$validated[$i+1];
	#Three different cases, first, internal, and last

	#First match (starts from start, does not look at previous
	if ($i==0) {
	    $ntstart = $match->{ntstart};
	    $ntend = $next->{ntstart};
	    $cl = $match->{clone};
	    if ($match->{strand} == '+1') {
		$clstart = $match->{clstart};
		$clend = $match->{clend} - ($match->{ntend} - $next->{ntstart});
	    }
	    else {
		$clstart = $match->{clstart} + ($match->{ntend} - $next->{ntstart});  
		$clend = $match->{clend};
	    }

	    $strand = $match->{strand};

	    if ($ntstart !=1 ) {
		print STDERR "BUG: ntstart of first match for $nt-$cl is not 1!\n";
	    }
	}

	#Last
	elsif ($i == ($size-1)) {
	    $ntstart = $match->{ntstart}+1;
	    $ntend = $match->{ntend};
	    $cl = $match->{clone};
	    
	    if ($match->{strand} == '+1') {
		$clstart = $match->{clstart}+1;
		$clend = $match->{clend};
	    }
	    else {
		$clstart = $match->{clstart};
		$clend = $match->{clend}-1;
	    }
	    $strand = $match->{strand};
	    
	    #Check if the last match extends until the end 
	    #(bug fix for ambiguity characters in the end)
	    
	    if ($ntseq->length > $ntend) {
		#print STDERR "Detected missing tail\n";
		my $diff = ($ntseq->length - $match->{ntend});
		print STDERR "Missing $diff bases from the end in $nt-$cl match, adding them\n";
		#Check if extension on clone feasible?
		$ntend+=$diff;
		$clend+=$diff;
	    }
	}

	#All others (internal)
	else {
	    $ntstart = $match->{ntstart}+1;
	    $ntend = $next->{ntstart};
	    $cl = $match->{clone};
	    if ($match->{strand} == '+1') {
		$clstart = $match->{clstart}+1;
		$clend = $match->{clend} - ($match->{ntend} - $next->{ntstart});

	    }
	    else {
		$clstart = $match->{clstart} + ($match->{ntend} - $next->{ntstart});  
		$clend = $match->{clend}-1;
	    }
	    $strand = $match->{strand};

	}
	print STDERR "$nt\t$hsp\t$ntstart\t$ntend\t$strand\tF\t$cl\t$clstart\t$clend\n";
	print GP "$nt\t$hsp\t$ntstart\t$ntend\t$strand\tF\t$cl\t$clstart\t$clend\n";
    }
}
close (GP);



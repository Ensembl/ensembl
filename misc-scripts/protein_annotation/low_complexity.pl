#! /usr/local/bin/perl
#
# Perl script make_other_rdb
#
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

make_other_rdb - DESCRIPTION 

This script populates the pfamseq table in the pfamdev relation datbase.

=head1 SYNOPSIS

make_other_rdb.pl <pfamseq fasta file>

Options:

-help : Show this help message
-files : Use the ultra-fast (but mysql specific) file import mechanism to load data


=head1 DESCRIPTION
    
This script populates the other_region table in the pfam relational database.
It should be run once at each major release (i.e. each time the underlying 
database changes).

Note: This script creates temporary files in the cwd

The following region types are currently supported:

low_complexity
coiled_coil
trans_membrane 
signal_peptide

				    
=head1 CONTACT
						    
Mail pfam@sanger.ac.uk with any queries


=cut
    

## Note: This code splits pfamseq into about 5 chunks
## before feeding each file to seg. This is because
## seg segmentation faults with large input files


# Let the code begin...



    
use strict;
use Getopt::Long;
use FileHandle;




my ( $help,
     $path,
     $use_files,
     $read_fh, 
     $write_fh,
     $fname,
     $total, 
     $rdb,
     $filecount,
     $entries,
     @pfamseqlist,
     @datfiles,
     @list);

&GetOptions( 'files=s' => \$use_files,
	     'help' => \$help);

my $pfamseqfasta = shift;

exec('perldoc', $0) if $help or not defined $pfamseqfasta;

$filecount = 1;


print STDERR "Splitting pfamseq...\n";
$filecount = 0;
print "open $pfamseqfasta \n";
open(_PFAMSEQ, $pfamseqfasta) or die "Fatal: Could not open $pfamseqfasta";
while(<_PFAMSEQ>) {
    /^\>/ and do {
	if (not ($entries % 35000)) {
	    close(TMP);
	    open(TMP, ">pfamseq.tmp$$.$filecount");
	    push @pfamseqlist, "pfamseq.tmp$$.$filecount";
	    $filecount++;
	}
	$entries++;
    };
    print TMP;
}
close(TMP);
#my ($pf_file) = $pfamseqfasta =~ /(chunk\.\d+)/;
# first, lets add the low_complexity information from seg

open(_SEGOUTPUT, ">/work1/birney/mongin/prot_annotation/outputs/low_compl.dat") || die "Can't open low_compl.dat";
close(_SEGOUTPUT);

open(_SEGOUTPUT, ">>/work1/birney/mongin/prot_annotation/outputs/low_compl.dat") || die "Can't open file low_compl.dat" ;
foreach my $pfamseq_tmp (@pfamseqlist) {
    print STDERR "Doing seg $pfamseq_tmp -l\n";
    @list = ();
    $total = 0;

    open (_SEG, "seg $pfamseq_tmp -l |") or die "Could not open the seg command - $!";
    while(<_SEG>) {
	my ($seqid, $from, $to, $score);
	
	/^>(\S+)\((\d+)\-(\d+)\)\s+complexity=(\S+)/ && do {
	    $seqid = $1;
	    $from = $2;
	    $to = $3;
	    $score = $4;
	    
	    if ($use_files) {
		if (not defined $write_fh) {
		    $write_fh = FileHandle->new();
		    $fname = "$use_files/other_table$$";
		    $write_fh->open(">$fname") or do {
			warn "Could not open $fname for writing so skipping";
			next;
		    };
		}
	my $blank = "BLANK";
		print $write_fh "$seqid\t$from\t$to\t15\t0\t0\tlow_complexity\t0\t0\t$score\n";

		#print $write_fh "\\N\t$seqid\t$from\t$to\tlow_complexity\tseg\t$score\n";
	    }
	    else {

#Print in the following format to allow an easy loading into Ensembl database
#FORMAT: translation\tseq_start\tseq_end\tanalysis\thstart\thend\thid\tscore\tevalue\tperc_id\n
			print _SEGOUTPUT "$seqid\t$from\t$to\t15\t0\t0\tlow_complexity\t0\t0\t$score\n";

			#print _SEGOUTPUT "\\N\t$seqid\t$from\t$to\tlow_complexity\tseg\t$score\n";


		push @list, { 'seqid' => $seqid,
			      'from' => $from,
			      'to' => $to,
			      'type' => "low_complexity",
			      'source' => "seg",
			      'score' => $score
				   };
	    }
	    
	};
    }
    
    
    eval {
	if ($use_files) {
	    $write_fh->close or die "Could not close $fname for load";

	    $write_fh = undef;
	}
	else {
	  
	}
    };
    $@ and die "Fatal error: Failed to add records to pfamseq table from $fname [$@]";
}
close(_SEGOUTPUT);


# clean up...

foreach my $tmpfile (@pfamseqlist) {
    unlink $tmpfile;
}
if ($use_files) {
    unlink $fname;
}

print STDERR "Total rows inserted = $total\n";
exit(0);



#!/usr/local/bin/perl

=head1 fpc2staticgp

This script turns a modified agp file (with fpc information)
into a tab-delimited file for the static golden path table.

you have to dump out a file of

  contig-id,internal_id,start,length 

from the mysql database you are loading into. (ie, go

  select id,internal_id,offset,length from contig into OUTFILE 'somewhere';

The script is run as

  perl fpc2staticgp dump.file *.fpc > static.txt

To make the fpc file, use mergechrfpc.pl

=cut


my $dump = shift;


open(L,">parse.log");
open(D,$dump) || die "Could not open $dump\n";
while(<D>) {
    ($id,$internal_id,$start,$end) = split;
    $id =~ /(\S+)\./ || die "Bad id $id";
    $clone = $1;
    if( !defined $h{$clone} ) {
	$h{$clone} = [];
    }

    push(@{$h{$clone}},$internal_id);
    $start{$internal_id} = $start;
    $end{$internal_id} = $start+$end-1;
}

print STDERR "Hash is loaded\n";

foreach my $f ( @ARGV ) {
    open(F,$f) || die "could not open $f\n";
    while( <F> ) {
	# pick up only contig lines. have a +/-
	/\-|\+/ || do { next; };
	my ($chr,$start,$end,$number,$d,$id,$rstart,$rend,$ori,$fpc,$fstart,$fend) = split;
	if( !defined $fend ) {
	    die "on $_\n";
	}

	($clone) = $id =~ /(\S+)\./;
	#print STDERR "Looking at $clone $id\n";
	$seen = 0;
	foreach $iid ( @{$h{$clone}} ) {
	    #print STDERR "Looking at $iid $rstart:$rend ",$start{$iid}," ",$end{$iid},"\n";

	    if( $rstart >= $start{$iid} && $rstart <= $end{$iid} ) {
		if( $rend > $end{$iid} ) {
		    print L "$f $clone ($id) Bad news... $iid does not fit! $rend vs ",$end{$iid},"\n";
		} else {
		    $seen = 1;
		    if( $ori eq '-' ) {
			$orit = -1;
		    } else {
			$orit = 1;
		    }
		    $rstart = $rstart - $start{$iid} +1;
		    $rend   = $rend   - $start{$iid} +1;


		    print "$fpc\t$chr\t$iid\t$start\t$end\t$fstart\t$fend\t$rstart\t$rend\t$orit\tUCSC\n";
		    last;
		}
	    }
	}
	if( $seen == 0 ) {
	    print STDERR "Unable to fit $clone $rstart:$rend\n";
	}
    }
}
		   
		   
		

    

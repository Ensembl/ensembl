my $old;
my $new;
while (<>) {
    
    if (/grep (\S+)/) {
	$old=$1;
	next;
    }

    if (/(ENSE\S+)\s(ENSE\S+)/) {
	if ($2 eq $old) {
	    $new=$1;
	    if ($old ne $1) {
		print "update translation set start_exon\=\"$1\" where start_exon\=\"$old\"\;\n";
		print "update translation set end_exon\=\"$1\" where end_exon\=\"$old\"\;\n";
		print "update supporting_feature set exon\=\"$1\" where exon\=\"$old\"\;\n";
		print "update exon_transcript set exon\=\"$1\" where exon\=\"$old\"\;\n";
	    }
	}
	elsif ($1 eq $old) {
	    if ($new ne $2) {
		print "update translation set start_exon\=\"$new\" where start_exon\=\"$2\"\;\n";
		print "update translation set end_exon\=\"$new\" where end_exon\=\"$2\"\;\n";
		print "update supporting_feature set exon\=\"$new\" where exon\=\"$2\"\;\n";
		print "update exon_transcript set exon\=\"$new\" where exon\=\"$2\"\;\n";
	    }
	    $old = $2;
	}
    }
}

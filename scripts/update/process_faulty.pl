while (<>) {
    if (/Exception in (\.*)/) {
	print "$1\n";
    }
}

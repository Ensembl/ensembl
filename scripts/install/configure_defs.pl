use strict;

foreach (@ARGV) {
	reconfigure($_);
}

sub reconfigure {
	my $file = shift;
	open IN, $file or die "Cant open $file: $!\n";
	local $/ = undef;
	my $txt = <IN>;
	close IN;
	foreach my $k (keys %ENV) {
	    $k =~ /^ENSEMBL_/ || next;
	    $txt =~ s/(\$${k}\s*?=).*?\n/$1 '$ENV{$k}';\n/gs;
	}
	open OUT, ">$file" or die "Cant write to $file: $!\n";
	print OUT $txt;
	close OUT;
}

my @files;
while (<>) {
	chomp();
	push (@files,$_);
}
FILE: foreach my $file (@files) {
	open (FILE,$file);
	if ($file =~ /val\=(NT_\d+)/) {
		$nt =$1;
	}
	while (<FILE>) {
	    if (/CONTIG/) {
		$line =$_;
		$line =~ s/\<a href\S+val\=\S+\.\d+\>//g;
		$line =~ s/\<\/a\>//g;
		$line =~ s/CONTIG//g;
		print "$nt assembly:\n";
		print $line;
		while (<FILE>) {
		    $line =$_;
		    $line =~ s/\<a href\S+val\=\S+\.\d+\>//g;
		    $line =~ s/\<\/a\>//g;
		    if ($line =~ /\<\/pre/) {next FILE;}
		    if ($line =~ /\/\//) {next FILE;}
		    print $line;
		}
	    }
	}
	print "\n";
	close(FILE);
}

my %map;
my $file = pop (@ARGV);
foreach my $f (@ARGV) {
    #print "Reading $f\n";
    open ($f,"<$f");
    while (<$f>) {
	chomp;
	$_ =~ /(\S+)\s+(\S+)/;
	$map{$2}=$1;
    }
    close ($f);
}
open (FILE,"<$file");
while (<FILE>) {
    chomp;
    my @tabs = split(/\t/);
    my $line;
    foreach my $tab (@tabs) {
	if ($map{$tab}) {
	    $line .= $map{$tab}."\t";
	}
	else {
	    $line .= $tab."\t";
	}
    }
    $line =~ s/\t$//g;
    print "$line\n";
}

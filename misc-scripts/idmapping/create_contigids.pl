my %contigs;
while (<>) {
    chomp;
    if ($_ =~ /(\S+)\.\d+\.\d+\.\d+/) {
	push (@{$contigs{$1}},$_);
    }
    elsif ($_ =~ /(\S+)\.\d+\.(\d+)/) {
	if ($2 =~ /^1/) {
	    push (@{$contigs{$1}},$_);
	}
	else {
	    #print "Couldn't parse line $_\n";
	} 
    }
    else {
	#print "Couldn't parse line $_\n";
    }
}
my %idmap;
foreach my $clone (keys %contigs) {
    my %clone;
    my $c=0;
    foreach my $line (@{$contigs{$clone}}) {
	#print "Parsing line $line\n";
	if ($line =~ /\S+\.\d+\.(\d+)\.\d+/) {
	    $clone{$1}=$line;
	}
	elsif ($_ =~ /(\S+)\.\d+\.(\d+)/) {
	    $clone{1}=$line;
	    $c++;
	}
	if ($c > 1) {
	    #print "BARF\n";
	}
    }
    #print "Sorted by start lines for clone $clone:\n";
    my $c=1;
    foreach $line (sort {$a <=> $b} (keys %clone)) {
	my $id = "0000$c";
	while (length ($id) > 5) {
	    $id =~ s/^0//g;
	}
	my $final = "$clone.$id";
	$idmap{$line}=$final;
	print $clone{$line}." $final\n";
	$c++;
    }
    #print "\n";
}
	    
	

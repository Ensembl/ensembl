use strict;

my ($input)= @ARGV;

open (IN,"$input") || die "Can't open input file: $input\n";

my %hash;

my %freq;

while (<IN>) {
    chomp;
    my @array = split;

    my $dna_id = $array[0];
    my $start = $array[1];
    
    push (@{$hash{$dna_id}},$array[1]);

}

foreach my $k(keys %hash) {
    #print "HERE\n";
    my $w_start = 1;
    my $w_end = 500;
    
    my @array = @{$hash{$k}};

    @array =  sort {$a <=> $b} @array;

    my $scalar = scalar @array;
    my $max = $array[$scalar - 1];
    
       
    while ($w_end < $max) {
	my $count;
	foreach my $start (@{$hash{$k}}) {
	    if (($start >= $w_start) && ($start <= $w_end)) {
		$count++;
		$freq{$count}++;
	    } 
	}
	$w_start = $w_start + 500;
	$w_end = $w_end + 500;
	
	if ($freq{$count}) {
	    #print  "$k\t$count\t$freq{$count}\n";
	}
    }
}

foreach my $a (keys %freq) {
    print "$a\t$freq{$a}\n";
}

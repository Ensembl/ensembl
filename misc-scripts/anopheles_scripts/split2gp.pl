

use strict;


open (GOLDEN,"/acari/work7a/mongin/assembly/golden_path_base.txt") || die;


my %celera_con;
my %celera_len;
my %seen;


while( <GOLDEN> ) {

    my($celera,$raw_id,$fstart,$fend,$contig_start,$contig_end) = split;

   #$celera =~ s/^CHR//g;
    my $h = {};

    $h->{'raw_id'} = $raw_id;
    $h->{'fstart'} = $fstart;
    $h->{'fend'} = $fend;
    $h->{'cstart'} = $contig_start;
    $h->{'cend'} = $contig_end;


    if( !defined $celera_con{$celera} ) {
	$celera_con{$celera} = [];
    }

    push(@{$celera_con{$celera}},$h);

    if( $celera_len{$celera} < $fend ) {
	$celera_len{$celera} = $fend;
    }
}


open (IN,"/acari/work7a/mongin/assembly/Agambiae_golden_path_17_7_02_XL_formated.txt");

while (<IN>) {

    my ($celera_con,$chr,$start,$end,$ori) = split;

    $seen{$celera_con} = 1;

    if( $ori eq '+' ) {
	foreach my $h ( @{$celera_con{$celera_con}} ) {
	    print $celera_con,"\t",$chr,"\t",$h->{'raw_id'},"\t",$start+$h->{'fstart'}-1,"\t",$start+$h->{'fend'}-1,"\t",$h->{'fstart'},"\t",$h->{'fend'},"\t1\t",$h->{'cstart'},"\t",$h->{'cend'},"\t1\tCELERA\n";
	}
    } else {
	foreach my $h ( reverse @{$celera_con{$celera_con}} ) {
	    my $len = $celera_len{$celera_con};
	    my $local_start = $len - $h->{'fend'} + $start;
	    my $local_end = $len - $h->{'fstart'} + $start;

	    #print $h->{'raw_id'},"\t",$local_start,"\t",$local_end,"\t-1\t",$h->{'cstart'},"\t",$h->{'cend'},"\n";
	    print $celera_con,"\t",$chr,"\t",$h->{'raw_id'},"\t",$local_start,"\t",$local_end,"\t",$h->{'fstart'},"\t",$h->{'fend'},"\t-1\t",$h->{'cstart'},"\t",$h->{'cend'},"\t-1\tCELERA\n";
	}
    }
}


# foreach celera contig, if not seen, then put into unknown.

my $start = 1;




foreach my $con ( keys %celera_con ) {
    if( $seen{$con} ) {
	next;
    }

    
    foreach my $h ( @{$celera_con{$con}} ) {
	print $con,"\t",'UNKN',"\t",$h->{'raw_id'},"\t",$start+$h->{'fstart'}-1,"\t",$start+$h->{'fend'}-1,"\t",$h->{'fstart'},"\t",$h->{'fend'},"\t1\t",$h->{'cstart'},"\t",$h->{'cend'},"\t1\tCELERA\n";
    }
    
    $start += $celera_len{$con} + 1000;
}




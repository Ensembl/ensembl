use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;

use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

my $host      = 'localhost';
my $dbuser    = 'manu';
my $dbname    = 'anopheles_gambiae_core_10_2';
my $dbpass    = '';
my $path      = 'MOZ2';

my %hash;
my $total;


print STDERR "Connecting to $host, $dbname\n";


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );

my $clone_adapt = $db->get_CloneAdaptor();
my $slice_adapt = $db->get_SliceAdaptor();

my $query1 = "select clone_id,name from clone";
my $sth1 = $db->prepare($query1);
$sth1->execute();

while(my ($id,$clone_name) = $sth1->fetchrow_array) {

    my $slice = $slice_adapt->fetch_by_clone_accession($clone_name);

    my $chr_name = $slice->chr_name;

    my $seq = $slice->seq;
    
    &get_single($seq);
    &get_pairs($seq);

}


print "A: ".$hash{'A'}/$total."\n";
print "C: ".$hash{'C'}/$total."\n";
print "T: ".$hash{'T'}/$total."\n";
print "G: ".$hash{'G'}/$total."\n";

print "AG: ".$hash{'AG'}/$total."\n";
print "AC: ".$hash{'AC'}/$total."\n";
print "AT: ".$hash{'AT'}/$total."\n";
print "AA: ".$hash{'AA'}/$total."\n";
print "CA: ".$hash{'CA'}/$total."\n";
print "CC: ".$hash{'CC'}/$total."\n";
print "CT: ".$hash{'CT'}/$total."\n";
print "CG: ".$hash{'CG'}/$total."\n";
print "TA: ".$hash{'TA'}/$total."\n";
print "TC: ".$hash{'TC'}/$total."\n";
print "TT: ".$hash{'TT'}/$total."\n";
print "TG: ".$hash{'TG'}/$total."\n";
print "GA: ".$hash{'GA'}/$total."\n";
print "GC: ".$hash{'GC'}/$total."\n";
print "GT: ".$hash{'GT'}/$total."\n";
print "GG: ".$hash{'GG'}/$total."\n";



sub get_single {
    my ($seq) = @_;
    
    my @seqarray = split('',$seq);

    foreach my $base (@seqarray) {
	if ($base =~ /A/) {
	    $hash{'A'}++;
	}
	if ($base =~ /C/) {
	    $hash{'C'}++;
	}
	if ($base =~ /T/) {
	    $hash{'T'}++;
	}
	if ($base =~ /G/) {
	    $hash{'G'}++;
	}
    }
}

sub get_pairs {
    my ($seq) = @_;
    my $ind = 1;
    my $length = length($seq);
    
    while($ind<=$length) {
	my $pair = substr($seq,$ind,2);
	
	if ($pair =~ /^A/i) {
	    if ($pair =~ /AA/i) {
		$hash{'AA'}++;
	    }
	    elsif ($pair =~ /AC/i) {
		$hash{'AC'}++;
	    }
	    elsif ($pair =~ /AT/i) {
		$hash{'AT'}++;
	    }
	    elsif ($pair =~ /AG/i) {
		$hash{'AG'}++;
	    }
	}

	elsif ($pair =~ /^C/i) {
	    if ($pair =~ /CA/i) {
		$hash{'CA'}++;
	    }
	    elsif ($pair =~ /CC/i) {
		$hash{'CC'}++;
	    }
	    elsif ($pair =~ /CT/i) {
		$hash{'CT'}++;
	    }
	    elsif ($pair =~ /CG/i) {
		$hash{'CG'}++;
	    }
	}
	
	elsif ($pair =~ /^T/i) {
	    if ($pair =~ /TA/i) {
		$hash{'TA'}++;
	    }
	    elsif ($pair =~ /TC/i) {
		$hash{'TC'}++;
	    }
	    elsif ($pair =~ /TT/i) {
		$hash{'TT'}++;
	    }
	    elsif ($pair =~ /TG/i) {
		$hash{'TG'}++;
	    }
	}
	
	elsif ($pair =~ /^G/i) {
	    if ($pair =~ /GA/i) {
		$hash{'GA'}++;
	    }
	    elsif ($pair =~ /GC/i) {
		$hash{'GC'}++;
	    }
	    elsif ($pair =~ /GT/i) {
		$hash{'GT'}++;
	    }
	    elsif ($pair =~ /GG/i) {
		$hash{'GG'}++;
	    }
	}
	$ind++;
    }
    
    $total = $total + $ind;
    
}

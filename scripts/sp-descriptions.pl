#!/usr/local/bin/perl
# $Header$

# Create a file dbname\taccno\tdescription from swissprot flat file.
# Used by gene-descriptions.pl

# Usage: sp-descriptions.pl < swiss-prot-flatfile  > sp-descriptions.dat

my $acc; 
my $desc="";
while( <> ) {
    if ( /^ID/ ) {
        if (/PRELIMINARY;/) {
            $db='SPTREMBL';
        } elsif (/STANDARD;/) { 
            $db='SWISS-PROT';
        } else {
            chomp;
            die "can't recognize: $_";
        }
    }
    if ( /^AC\s+(\S+)/ )  {
        $acc = $1;
        $acc =~ s/;$//g;
    }

    if ( /^DE\s+(\S.*\S)/) {
        $desc .= $1;
    }

    if ( m://: ) {
        $desc =~ s/\{.*\}//g; 
        $sp_desc{"$db\t$acc"} = $desc; 
        $acc = undef; 
        $desc = ""; 
    }
}

foreach $acc ( keys %sp_desc )  {
    $desc = $sp_desc{$acc};
    $desc =~ s/\n/\\n/g;
    print "$acc\t$desc\n";
}

exit;

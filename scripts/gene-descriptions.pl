#!/usr/local/bin/perl
# $Header$
# 
# Script to create a tab-separated file (to be loaded into the
# gene_description table) from a file sp-descriptions.dat
# (swiss-prot/sptrembl -> description mapping, from running
# sp-descriptions.pl on a swissprot/trembl flatfile) and from the
# mapping.dat (ensg -> sp/trembl accno mapping, obtained from script
# sp-mapping-dump.sh). 


## list of regexps tested, in order of increasing preference, when
## deciding which description to use (used by sub compare_desc):
my @word_order = 
   qw(unknown hypothetical putative novel probable [0-9]{3} kDa fragment cdna protein);

$Usage = "Usage: $0 sp-descriptions.dat mapping.dat > gene-descriptions.dat\n";

die $Usage if @ARGV == 0;

$spdesc = shift;
$map = shift;

open(SPDESC, $spdesc) || die "$spdesc:$!";

undef %sp_desc;
while(<SPDESC>) {
    chomp;
    ($db, $acc, $desc)=split(/\t/);

    $sp_desc{"$db:$acc"}=$desc;
}

close(SPDESC) || die "$!";

open(MAP, $map) || die "$map:$!";

undef %gene_desc;

LINE:
while (<MAP>) {                       
    chomp;
    ($ensp, $ensg, $db, $acc)=split(/\t/);

    if ( defined($gene_desc{$ensg}) ) {
        ($prevdb, $prev_desc)  = @{$gene_desc{$ensg}};
        if ($prevdb  eq 'SWISS-PROT') {   
            next LINE;                  # nothing to change
        }
        if ($db  eq 'SWISS-PROT') {   
            $desc = $sp_desc{"$db:$acc"};
            $gene_desc{$ensg} = [ $db, $desc]; # kick out the SPTREMBL desc.
            next LINE;
        }

        if ($db  eq 'SPTREMBL' &&  $prevdb eq $db ) {   
            $desc = $sp_desc{"$db:$acc"}; 
            if ( &compare_desc($prev_desc, $desc) < 0 ) {
                # new desc is better
                # warn "new better: $desc (old was: $prev_desc)\n";
                $gene_desc{$ensg} = [ $db, $desc];
                next LINE;
            } else {
                # warn "old better: $prev_desc (new is: $desc)\n";
                next LINE;
            }
            die "should not reach this point: only know SWISS-PROT and SPTREMBL";
        }
    } else {
        $desc = $sp_desc{"$db:$acc"};
        $gene_desc{$ensg} = [ $db, $desc];
    }
}                                       # while <MAP>

#  now dump the stuff to stdout.
foreach $ensg ( keys %gene_desc )  { 
    ($db, $desc)   = @{$gene_desc{$ensg}};
    print STDOUT "$ensg\t$desc\n";
}

#### following taken from ensembl-external/scripts/family-input.pl

sub compare_desc {
    my ($a, $b) = @_; 

    my ($am, $bm);
    foreach my $w (@word_order) {
        $am = ($a =~ /$w/i)?1:0;
        $bm = ($b =~ /$w/i)?1:0;

        if ($am  != $bm ) {             # ie, one matches, other doesn't
            if ( $am == 1 ) {           # first one worse than second 
                return -1;
            } else { 
                return 1; 
            }
        }
    }
    # still look same; base result on length: longer is better
    return length($a) <=> length($b);
}                                       # compare_desc

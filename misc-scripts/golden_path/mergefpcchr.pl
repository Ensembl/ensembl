#!/usr/local/bin/perl


=head1 mergefpcchr

This script merges FPC data with chromosomal data agps. It relies heavily
on jim kents directory structure which looks like

For chromosomes

    chr-directory/chr_XXX_YYY.agp
  
For contigs

    contig-directory/chr_XXX/ctgZZZ/ctgZZZ.agp

script goes 

     perl mergefpcchr.pl chr-directory contig-directory

and builds a new chr_XXX_YYY.agp.fpc with the fpc information merged into the
agp information

=cut


my $chr  = shift;
my $contigdir = shift;

opendir(CHR,$chr) || die "no chromosome directory";
@chromo = readdir(CHR);
closedir(CHR);
shift(@chromo);
shift(@chromo);

foreach $chromo ( @chromo ) {
    #$chromo =~ /.*?(\d+).*\.agp/ || next;
    $chromo =~ /chr([^_\.]+)[_|\.].*agp/ || next;
    $chrnumber = $1;

    print STDERR "Doing $chromo\n";
    #if ($chromo ne 'chr22.agp') {
    #	next;
    #}
    open(GF,"$chr/$chromo") || die "no $file $!\n";
    open(OUT,">$chr/$chromo.fpc");
    

    opendir(DIR,"$contigdir/$chrnumber") || die "no $contigdir/$chrnumber $!\n";
    @contigs = readdir(DIR);
    closedir(DIR);
    
    shift(@contigs);
    shift(@contigs);
    

    foreach $contig ( @contigs ) {
	open(F,"$contigdir/$chrnumber/$contig/$contig.agp") || die "did not open $contig $!\n";
	while(<F>) {
	    # chr20   1       1970    1       P       AL360078.3      76699   78668   -
	    /\S+\/(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+\S\s+(\S+)\s+(\d+)\s+(\d+)/ ||next;
	    
	    $id = $1;
	    $fstart = $2;
	    $fend = $3;
	    $accf = $4;
	    $accst = $5;
	    $accend = $6;
	    $idstring = "$accf:$accst:$accend";
	    
	    #print "Storing contig string $idstring\n";
	    $fpc{$idstring} = $id;
	    $fpcstart{$idstring} = $fstart;
	    $fpcend{$idstring} = $fend;
	}
    }
    
    while( <GF> ) {
	#19/ctg113	1	37401	1	F	AC011523.3	1	37401	+
	#22/chr22       1       37693   1       F       AP000522.1      1       37693   +
	/(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+\S\s+(\S+)\s+(\d+)\s+(\d+)/ || do { print OUT $_; next; };
	
	$id = $1;
	$fstart = $2;
	$fend = $3;
	$accf = $4;
	$accst = $5;
	$accend = $6;
	$idstring = "$accf:$accst:$accend";
	#print "Reading chromosome line with $idstring\n";
	if( !defined $fpc{$idstring} ) {
	    print STDERR "Could not find fpc contig for $idstring\n";
	    next;
	}
	chomp;
	print OUT "$_ $fpc{$idstring} $fpcstart{$idstring} $fpcend{$idstring}\n";
    }
}


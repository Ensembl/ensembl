use strict;
use Getopt::Long;

my $file;

&GetOptions(
            'file:s'=>\$file
	    );


&split_file();
&run_low_complexity();
&run_coils();
&run_tmhmm();
&run_sign_pep();

sub split_file {
#Split Ensembl fasta file into chiunks of 2000 sequences
print STDERR "Splittinq $file\n";
open (ENSEMBL, "$file");
my $count = 0;
my $chunk = 1;

    $/ = "\>";
    
    while(<ENSEMBL>){
	
	if ($_ ne "\>") {
	    if ($count == 0) {
		open (CHUNK,">chunks/chunk.$chunk");
	    }
	    
	    $_ =~ s/\>$//;  
	    
	    print CHUNK ">$_";
	    $count++;
	    if ($count == 1000) {
		$count = 0;
		$chunk++;
	    }
	}
    }
}

sub run_low_complexity {
    print STDERR "Sending low_complexity job to the LSF queue\n";
    my $run_low_compl = "bsub -q acari -o tmp/low_compl.out -e tmp/low_compl.err perl /work1/birney/mongin/src/ensembl-live/misc-scripts/protein_annotation/low_complexity.pl $file";
    
    system($run_low_compl)==0 || die "$0\Error running '$run_low_compl' : $!";
}


sub run_coils {
    print STDERR "Sending coils jobs  to the LSF queue\n";
    my $dir = "chunks/";
    opendir(DIR,$dir) || die "No directory:$!";
    
    my @allfiles = readdir DIR;
    closedir DIR;
    
    foreach my $fasta_chunk(@allfiles) {
	if (($fasta_chunk ne ".") && ($fasta_chunk ne "..")) {
	    my $run_coils = "bsub -q acari -o tmp/$fasta_chunk.coils.out -e tmp/$fasta_chunk.coils.err perl /work1/birney/mongin/src/ensembl-live/misc-scripts/protein_annotation/coils.pl /work1/birney/mongin/prot_annotation/chunks/$fasta_chunk";
	    
	    system($run_coils)==0 || die "$0\Error running '$run_coils' : $!";
	}
    }
}

sub run_tmhmm {
    print STDERR "Sending tmhmm jobs  to the LSF queue\n";
    
    my $dir = "chunks/";
    opendir(DIR,$dir) || die "No directory:$!";
    
    my @allfiles = readdir DIR;
    closedir DIR;
    
    foreach my $fasta_chunk(@allfiles) {
	if (($fasta_chunk ne ".") && ($fasta_chunk ne "..")) {
	    my $run_tmhmm = "bsub -q acari -o tmp/$fasta_chunk.tmhmm.out -e tmp/$fasta_chunk.tmhmm.err perl /work1/birney/mongin/src/ensembl-live/misc-scripts/protein_annotation/tmhmm.pl /work1/birney/mongin/prot_annotation/chunks/$fasta_chunk";
	    
	    system($run_tmhmm)==0 || die "$0\Error running '$run_tmhmm' : $!";
	}
    }
}

sub run_sign_pep {
    print STDERR "Sending signal peptides jobs to the LSF queue\n";
    my @runfiles;
    my $dir = "chunks/";
    opendir(DIR,$dir) || die "No directory:$!";
    my @allfiles = readdir DIR;
    closedir DIR;
    
    my $sigp = "perl /work1/birney/mongin/src/ensembl-live/misc-scripts/protein_annotation/sigp.pl $file";
	    
    system($sigp)==0 || die "$0\Error running '$sigp' : $!";

    foreach my $file (@allfiles) {
	if ($file =~ /sigp_split.(\d+)/) {
	    my ($run) = $file =~ /(sigp_split.\d+)/;
	    print STDERR "RUN1: $run\n";
	    push(@runfiles,$file);
	}
    }

    foreach my $run (@runfiles) {
	print STDERR "RUN2: $run\n";
	my $run_peps = "bsub -q acari -o tmp/$run.peps.out -e tmp/$run.peps.err perl /work1/birney/mongin/src/ensembl-live/misc-scripts/protein_annotation/sigp_readsplit.pl $run";
	    
	 system($run_peps)==0 || die "$0\Error running '$run_peps' : $!";
     }
}

      		

   	    












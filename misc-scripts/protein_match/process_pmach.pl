use strict;

=head1 Process pmatch

=head1 Description

=head2 Aims
This script aims to run pmatch and postprocess pmatch to map Ensembl peptides to external databases (currently Swissprot and Refseq but may be extented). The first part of the script runs pmatch, the second part gets the percentage of a match of a unique Ensembl peptide which match to an unique external protein. The third part classify each ensembl match as PRIMARY match (the longest one and the one which will be used for the mapping, PSEUDO, DUPLICATE and REPEAT (pretty arbitrary criterias but may be useful for quality control).
NB: All of the intermediary files are written.

=head2 Options
-ens : Ensembl peptide fasta file
-sp : SP, SPTREMBL fasta file
-refseq: Refseq peptide fasta file

=head2 Contacts
mongin@ebi.ac.uk
birney@ebi.ac.uk

=cut  

use Getopt::Long;

my ($ens,$sp,$refseq,$pdb);

&GetOptions(
            
            'ens:s'=>\$ens,
            'sp:s'=>\$sp,
            'refseq:s'=>\$refseq,
	    'pdb:s'=>\$pdb
            );

&runpmatch();
&postprocesspmatch($sp);
&postprocesspmatch($refseq);
<<<<<<< process_pmach.pl

=======
&postprocesspmatch($pdb);
>>>>>>> 1.8
&finalprocess($sp);
&finalprocess($refseq);
&finalprocess($pdb);

#perl ../../../src/ensembl-live/misc-scripts/protein_match/process_pmach.pl -ens ../primary/SPAN_pepfile  -sp ../primary/SPTr.human.expanded -refseq ../primary/hs2.fsa -pdb ../primary/scop_human.fas 

sub runpmatch {
    print STDERR "Running pmatch\n";
 
#Run pmatch and store the data in files which will be kept for debugging
    my $pmatch1 = "/nfs/griffin2/rd/bin.ALPHA/pmatch -T 14 $sp $ens > ens_sp_rawpmatch";
    my $pmatch2 = "/nfs/griffin2/rd/bin.ALPHA/pmatch -T 14 $refseq $ens > ens_refseq_rawpmatch";
    #my $pmatch3 = "/nfs/griffin2/rd/bin.ALPHA/pmatch -T 14 $pdb $ens > ens_pdb_rawpmatch";


    system($pmatch1); # == 0 or die "$0\Error running '$pmatch1' : $!";
    system($pmatch2); #== 0 or die "$0\Error running '$pmatch2' : $!";
    #system($pmatch3); #== 0 or die "$0\Error running '$pmatch2' : $!"; 

}


    
sub postprocesspmatch {
    my ($db) = @_;
    my %hash1;
    my %hashlength;
      
#Post process the raw data from pmatch
    if ($db eq $sp) {
	print STDERR "Postprocessing pmatch for SP mapping\n";
	open (OUT, ">ens_sp.processed") || die "Can't open File\n";
	open (PROC, "ens_sp_rawpmatch") || die "Can't open File\n";
    }
    
    elsif ($db eq $refseq) {
	print STDERR "Postprocessing pmatch for REFSEQ mapping\n"; 
	open (OUT, ">ens_refseq.processed") || die "Can't open File\n";;
	open (PROC, "ens_refseq_rawpmatch") || die "Can't open file ens_refseq_rawpmatch\n";
    }

    elsif ($db eq $pdb) {
	print STDERR "Postprocessing pmatch for PDB mapping\n"; 
	open (OUT, ">ens_pdb.processed") || die "Can't open File\n";;
	open (PROC, "ens_pdb_rawpmatch") || die "Can't open file ens_pdb_rawpmatch\n";
    }
    
    while (<PROC>) {
#538     COBP00000033978 1       538     35.3    Q14146  1       538     35.3 
	my ($len,$id,$start,$end,$tperc,$query,$qst,$qend,$perc) = split;

	if ($db eq $refseq) {
	    #Get only the refseq ac (NP_\d+) 
	    ($query) = $query =~ /\w+\|\d+\|\w+\|(\w+)/;
	}

	my $uniq = "$id:$query";
	
#Add the percentage of similarity for the Ensembl peptide for a single match
#There is a bug at this step, some similarities can be over 100% !!! This problem may be solved by changing pmatch source code
	$hash1{$uniq} += $perc;
	$hashlength{$uniq} += $len;
	
    }

#Write out the processed data
    foreach my $key ( keys %hash1 ) {
	#if (($hashlength{$key} >= 20)) {
	if (($hash1{$key} >= 25)) {
	    ($a,$b) = split(/:/,$key);
	    print OUT "$a\t$b\t$hash1{$key}\n";
	}
	#else {
	 #   print "$a\t$b\t$hash1{$key}\t$hashlength{$key}\n";
	#}
    }
    close (PROC);
    close (OUT);                               
}

sub finalprocess {
#This final subroutine will use the postprocessed pmatch file and get back the best Ensembl match (labelled as PRIMARY) for a given external known protein.
    my ($db) = @_;
    
    if ($db eq $sp) {
	print STDERR "Getting final mapping for SP mapping\n";
	open (PROC, "ens_sp.processed");
	open (OUT, ">ens_sp.final");
    }

    elsif ($db eq $refseq) {
	print STDERR "Getting final mapping for REFSEQ mapping\n";
	open (PROC, "ens_refseq.processed") || die "Can' open file ens_refseq.processed\n";
	open (OUT, ">ens_refseq.final");
    }

    elsif ($db eq $pdb) {
	print STDERR "Getting final mapping for PDB mapping\n";
	open (PROC, "ens_pdb.processed") || die "Can' open file ens_refseq.processed\n";
	open (OUT, ">ens_pdb.final");
    }


    my %hash2;
    while (<PROC>) {
	my ($ens,$known,$perc) = split;
	#if ($perc > 100) {
	 #   print "$ens\t$known\t$perc\n";
	#}
	

	if( !defined $hash2{$known} ) {
	    $hash2{$known} = [];
	}
	
#Each single external protein correspond to an array of objects dealing with the name and the percentage of similarity of the Ensembl peptide matching with the the known external protein.
	my $p= NamePerc->new;
	$p->name($ens);
	$p->perc($perc);
	
	push(@{$hash2{$known}},$p);
    }

    foreach my $know ( keys %hash2 ) {
	my @array = @{$hash2{$know}};
	@array = sort { $b->perc <=> $a->perc } @array;
	
#The Ensembl match to the known protein is labelled as PRIMARY and will be used later for the mapping 
	my $top = shift @array;

	#if ($top->perc >= 20) {
	    
	    print  OUT "$know\t",$top->name,"\t",$top->perc,"\tPRIMARY\n";
	    
	    foreach $ens ( @array ) {
		
		if( $ens->perc > $top->perc ) {
		    die "Not good....";
		}
	    }

#If there is more than 20 Ensembl peptides matching a single known protein, these Ensembl peptides are labelled as REPEAT 	
	    if (scalar(@array) >= 20) {
		foreach my $repeat (@array) {
		    if( $repeat->perc+1 >= $top->perc ) {
			print OUT "$know\t",$repeat->name,"\t",$repeat->perc,"\tDUPLICATE\n";
		    }
		    else {
			
			print OUT "$know\t",$repeat->name,"\t",$repeat->perc,"\tREPEAT\n";
		    }
		}
	    }
	
#If less than 20, either duplicate if percentage of identity close to the PRIMARY labelled as DUPLICATE or labelled as PSEUDO. DUPLICATEs can also be used for the mapping 
	    if (scalar(@array) < 20) {
		foreach my $duplicate (@array) {
		    if( $duplicate->perc+1 >= $top->perc ) {
			print OUT "$know\t",$duplicate->name,"\t",$duplicate->perc,"\tDUPLICATE\n";
		    }
		    else {
			print OUT "$know\t",$duplicate->name,"\t",$duplicate->perc,"\tPSEUDO\n";
		    }
		} 
	    }              
	} 
    #}
    close (PROC);
    close (OUT);
}

#Set of objects to deal with the script

package NamePerc;


sub new {
    my $class= shift;
    my $self = {};
    bless $self,$class;
    return $self;
}


=head2 name

 Title   : name
 Usage   : $obj->name($newval)
 Function:
 Returns : value of name
 Args    : newvalue (optional)


=cut

sub name{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'name'} = $value;
    }
    return $obj->{'name'};

}

=head2 perc

 Title   : perc
 Usage   : $obj->perc($newval)
 Function:
 Returns : value of perc
 Args    : newvalue (optional)


=cut

sub perc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'perc'} = $value;
    }
    return $obj->{'perc'};

}                    



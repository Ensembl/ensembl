#!/usr/local/bin/perl

BEGIN{
 unshift(@INC, '../bioperl-live');
 unshift(@INC, '../ensembl');

}
use strict;


#use EnsWeb;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::AceDB::Obj;
use Getopt::Long;
use Time::Local;

my $dbtype = 'rdb';
my $host   = 'localhost';
my $dbname = 'ensembl';
my $dbuser = 'root';
my $mask   = 0;
my $dbpass = undef;
my $dbuser = 'ensro';
my $port   = 3306;
my $all    = 0;

my $infile  = '';

$| = 1;

&GetOptions( 'dbtype:s'  => \$dbtype,
	     'host:s'    => \$host,
	     'port:n'    => \$port,
	     'dbname:s'  => \$dbname,
	     'dbpass:s'  => \$dbpass,
             'dbuser:s'  => \$dbuser,
	     'infile=s'  => \$infile,
	     'all'       => \$all,
	     );

my ($db1)    = get_db_handles();
my ($clone1) = get_cloneids($infile,@ARGV);

if ($all)  {
  my @ids = $db1->get_all_Clone_id;
  $clone1 = \@ids;
}

my $i = 0;

my ($dnafh,$maskedfh,$genscanfh) = make_filehandles($dbname);

for ($i = 0; $i <= $#$clone1; $i++) {
    print(STDERR "Processing clone " . $clone1->[$i] . "\n");
    
    eval {
	my $cl1    = $db1->get_Clone($clone1->[$i]);
	
	foreach my $contig ($cl1->get_all_Contigs) {
	    
	    get_genscan_peptides($contig,$genscanfh);
	    get_dna             ($contig,$dnafh,$maskedfh);
	}
    };
    if ($@) {
	print("Error processing clone: $@\n");
    }
}

close($dnafh);
close($maskedfh);
close($genscanfh);


##################################################################
#
# Subroutines
#
##################################################################

sub make_filehandles {
    my ($dbname) = @_;


    my $dnafile     = $dbname . ".dna.fa";
    my $maskedfile  = $dbname . ".masked.fa";
    my $genscanfile = $dbname . ".genscan.fa";

    open (DNA,">$dnafile");
    open (MAS,">$maskedfile");
    open (GEN,">$genscanfile");


    return (\*DNA,\*MAS,\*GEN);
}

sub get_dna {
    my ($contig,$dnafh,$maskedfh) = @_;

    my $dnastr  = $contig->primary_seq->seq;
    my @repeats = $contig->get_all_RepeatFeatures;
	
    my $maskeddnastr  = mask_features($dnastr,@repeats);

    $dnastr        =~ s/(.{72})/$1\n/g;
    $maskeddnastr  =~ s/(.{72})/$1\n/g;
	    
    print($dnafh ">" .$contig->id . "\n$dnastr\n");
    print($maskedfh ">" .$contig->id . "\n$maskeddnastr\n");
}

sub get_genscan_peptides {
    my ($contig,$fh) = @_;

    my @genscan = $contig->get_all_PredictionFeatures();
     
    my $trancount = 1;
      
    foreach my $genscan (@genscan) {
	print STDERR "Got a genscan with ".$genscan->id."\n";


	my $transcript = new Bio::EnsEMBL::Transcript;
	$transcript->id($contig->id . "." . $genscan->seqname);
	$trancount++;

	my @exons;
	my $count    = 1;

	foreach my $f ($genscan->sub_SeqFeature) {
	    #print STDERR "A sub seq feature found..\n";
	    my $exon  = new Bio::EnsEMBL::Exon;
	    $exon->id       ($contig->id . ".$count");
	    $exon->contig_id($contig->id);
	    $exon->start    ($f->start);
	    $exon->end      ($f->end  );
	    $exon->strand   ($f->strand);
	    $exon->attach_seq($contig->primary_seq);
	    
	    push(@exons,$exon);
	    $count++;
	    
	}
	
	my $translation = new Bio::EnsEMBL::Translation;
	$translation->id($contig->id.".".$genscan->seqname);
	
	if ($exons[0]->strand == 1) {
	    @exons = sort {$a->start <=> $b->start} @exons;
	    $translation->start        ($exons[0]->start);
	    $translation->end          ($exons[$#exons]->end);
	    
	} else {
	    @exons = sort {$b->start <=> $a->start} @exons;
	    $translation->start        ($exons[0]->end);
	    $translation->end          ($exons[$#exons]->start);
	    
	}
	
	$translation->start_exon_id($exons[0]->id);
	$translation->end_exon_id  ($exons[$#exons]->id);
	
	my $endphase = 0;
	
	foreach my $exon (@exons) {
	    
	  $exon->phase         ($endphase);
	  $transcript->add_Exon($exon);
	  $endphase = $exon->end_phase;

      }

	$transcript->translation($translation);
	
	my $seq = $transcript->dna_seq;
	
	my $found;
	
	my $seq0    = $seq->translate('*','X',0);
	my $seqstr0 = $seq0->seq; chop($seqstr0);

	if ($seqstr0 !~ /\*/) {
	    $found = $seqstr0;
	}

	my $seq1    = $seq->translate('*','X',2);
	my $seqstr1 = $seq1->seq; chop($seqstr1);

	if ($seqstr1 !~ /\*/) {
	  $found = $seqstr1;
	}

	my $seq2    = $seq->translate('*','X',1);
	my $seqstr2 = $seq2->seq; chop($seqstr2);

	if ($seqstr2 !~ /\*/) {
	  $found = $seqstr2;
	}

	if (defined($found) && (length($found) > 0)) {
	    $found =~ s/(.{72})/$1\n/g;
	    print($fh ">" .$transcript->id . "\n$found\n");
	} else {
	    $transcript->warn("Couldn't translate " . $transcript->id . "\n");
	}
    }
}

sub mask_features {
    my ($dnastr,@repeats) = @_;

    my $dnalen = length($dnastr);

    REP:foreach my $f (@repeats) {

	my $start  = $f->start;
	my $end    = $f->end;
	my $length = ($end - $start) + 1;


	if ($start < 0 || $start > $dnalen || $end < 0 || $end > $dnalen) {
	    print STDERR "Eeek! Coordinate mismatch - $start or $end not within $dnalen\n";
	    next REP;
	}

	$start--;

	my $padstr = 'N' x $length;

	substr ($dnastr,$start,$length) = $padstr;

    }

    return $dnastr;
}
	
sub get_cloneids {
    my ($infile,@ARGV) = @_;
    
    my @clone1;
    my @acc;
    
    if (defined($infile)) {
	open(IN,"<$infile");
	while (<IN>) {
	    chomp;
            my ($clone,$acc) = split (' ',$_);	
	    push(@clone1,$clone);
            push(@acc,$acc); 
	}
	close(IN);
    }

    while ($#ARGV >= 0) {
	my $clone1 = shift @ARGV;
	push(@clone1,$clone1);
    }
    return (\@clone1);
}



sub get_db_handles {
    my ($db1);

    if( $dbtype =~ 'ace' ) {
	$db1 = Bio::EnsEMBL::AceDB::Obj->new( -host => $host, 
					      -port => $port);
    } elsif ( $dbtype =~ 'rdb' ) {

	my $locator = "Bio::EnsEMBL::DBSQL::Obj/host=$host;" .
  	              "port=$port;"      .
		      "dbname=$dbname;"  .
		      "user=$dbuser;"    .
		      "pass=$dbpass";

	$db1 = Bio::EnsEMBL::DBLoader->new($locator);
	
    } else {
	die("$dbtype is not a good type (should be ace or rdb)");

    }

    return ($db1);
}

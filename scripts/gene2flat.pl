#!/usr/local/bin/perl

=head1 NAME

gene2flat

=head1 SYNOPSIS
 
  gene2flat ENSG00000012

=head1 DESCRIPTION

gene2flat produces a number of flat file outputs of the genes,
in particular the protein translation

=head1 OPTIONS


    -host    host name for database (gets put as host= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::DBAdaptor)

    -format    [pep/dump/transcript/exon] dump in peptides/info/dna/exon format

    -getall    all genes from the database [no applicable to timdb]

    -usefile   read in on stdin a list of gene ids, one gene id per line

    -verbose   print to STDERR on each gene to dump

    -help      displays this documentation with PERLDOC

    -chunk     chunk size to array (only use if know what this means)

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBLoader;
use Bio::SeqIO;

use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

use Getopt::Long;

my $dbtype = 'rdb';
my $host   = 'ecs1c';
my $port   = '';
my $dbname = 'idmap_apr01';
my $dbuser = 'ensadmin';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';

my $format  = 'transcript';
my $usefile = 0;
my $getall  = 0;
my $verbose = 0;
my $webdir = undef;
my $logerror = undef;
my $help =0;
my $chunk = 1;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'chunk:i'    => \$chunk,
	     'usefile'    => \$usefile,
	     'format:s'   => \$format,
	     'getall'     => \$getall,
	     'verbose'    => \$verbose,
	     'webdir:s'     => \$webdir,
	     'logerror:s' => \$logerror,
	     'h|help'     => \$help
	     );
my $db;

if ($help) {
    exec('perldoc', $0);
}
open (ERR,">/nfs/acari/elia/panic/log");
my $locator = "$module/host=$host;port=;dbname=$dbname;user=$dbuser;pass=$dbpass";
$db =  Bio::EnsEMBL::DBLoader->new($locator);

my @gene_id;

if( $usefile ) {
    while( <> ) {
	my ($g) = split;
	push(@gene_id,$g);
    }
} elsif ( $getall == 1 ) {
    @gene_id = $db->gene_Obj->get_all_Gene_id();
} else {
    @gene_id = @ARGV;
}

my $seqio;

if( $format eq 'pep' || $format eq 'transcript' || $format eq 'exon') {
    $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;
}

if( $format eq 'id' ) {
    foreach my $id ( @gene_id ) {
	print "$id\n";
    }
    exit(0);
}

if ($format eq 'webdump') {
  if (! -e $webdir) {
    mkdir($webdir, 0777) or die "Can't create '$webdir' : $!";
  }
}

$db->DESTROY;
$db=undef;

while ( @gene_id > 0 ) {
    my @chunk_list = splice(@gene_id,0,$chunk);

    if( $verbose ) {
	print ERR "Fetching @chunk_list\n";
    }

    eval {
	    my $locator = "$module/host=$host;port=;dbname=$dbname;user=$dbuser;pass=$dbpass";
	    $db =  Bio::EnsEMBL::DBLoader->new($locator);
            $db->static_golden_path_type('UCSC');
            my $sa = $db->get_StaticGoldenPathAdaptor();

	my @genes = $db->gene_Obj->get_array_supporting('none',@chunk_list);
	foreach my $gene ( @genes ) {
	    my $gene_id = $gene->id();
	    if( $format eq 'pep' ) {
		foreach my $trans ( $gene->each_Transcript ) {
		    # get out first exon. Tag it to clone and gene on this basis
		    my @exon = $trans->each_Exon;
		    my $fe = $exon[0];

                    my ($chr,$gene_start,$cdna_start) = find_trans_start($trans);

		    my $tseq = $trans->translate();
		    if ( $tseq->seq =~ /\*/ ) {
			print ERR "translation of ".$trans->id." has stop codons. Skipping! (in clone". $fe->clone_id .")\n";
			#next;
		    }
                    my $gene_version = $gene->version;

		    $tseq->desc("Gene:$gene_id.$gene_version Clone:".$fe->clone_id . " Contig:" . $fe->contig_id . " Chr: " . $chr . " Pos: " . $cdna_start);
		    $seqio->write_seq($tseq);
		}
	    } elsif ( $format eq 'dump' ) {
		foreach my $trans ( $gene->each_Transcript ) {
		    print "Transcript ",$trans->id,"\n";
		    foreach my $exon ( $trans->each_Exon ) {
			print "  Exon ",$exon->id," ",$exon->contig_id,":",$exon->start,"-",$exon->end,".",$exon->strand,"\n";
			my $seq = $exon->seq();
			my $str = $seq->seq();
			print "    Start phase ",$exon->phase,"[",substr($str,0,10),"] End phase ",$exon->end_phase," [",substr($str,-10),"]\n";
		    }
		}
		
	    } 
	    elsif ($format eq 'exon') {
		foreach my $trans ( $gene->each_Transcript ) {
		    foreach my $exon ($trans->each_Exon) {
			my $seq = $exon->translate();
			$seq->id($exon->id);
			$seqio->write_seq($seq);
		    }
		}
	    }
	    elsif ($format eq 'transcript') {
		foreach my $trans ( $gene->each_Transcript ) {
                    my ($chr,$gene_start,$cdna_start) = find_trans_start($trans);
		    my $seq = $trans->dna_seq();
		    $seq->id($trans->id);
		    my @exon = $trans->each_Exon;
		    print ERR "Translation start".$trans->translation->start."\n";
		    my $fe = $exon[0];
                    my $gene_version = $gene->version;
		    $seq->desc("Gene:$gene_id.$gene_version Clone:".$fe->clone_id . " Contig:" . $fe->contig_id . " Chr: " . $chr . " Pos: " . $cdna_start);
		    $seqio->write_seq($seq);
		}
	    }
	    elsif ($format eq 'webdump') {
	      TRAN: foreach my $trans ( $gene->each_Transcript ) {

		  my $pep_file =  $webdir.$trans->id.".pep";

		  if (-e $pep_file) {
		    print ERR "Peptide file exists - skipping\n";
		    next TRAN;
		  }

		    my $trans_file = $webdir.$trans->id.".trans";
		    open (TRANS,">$trans_file");
		    my $seqiot = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*TRANS );
		    my $seq = $trans->dna_seq();
		    print ERR "dumping transcript",$trans->id,"\n";
		    $seq->display_id($trans->id);
		    my @exon = $trans->each_Exon;
		    my $fe = $exon[0];

                    my ($chr,$gene_start,$cdna_start) = find_trans_start($trans);
                    my $gene_version = $gene->version;
		    $seq->desc("Gene:$gene_id.$gene_version Clone:".$fe->clone_id . " Contig:" . $fe->contig_id . " Chr: " . $chr . " Pos: " . $cdna_start);
		    $seqiot->write_seq($seq);
		    $seqiot=undef;
		    close (TRANS);
		    

		    open (PEP,">$pep_file");
		    my $seqiop = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*PEP) ;
		    my $tseq = $trans->translate();
		    if ( $tseq->seq =~ /\*/ ) {
			print ERR "Skipping peptide dumping of ".$gene->id.", translation has stop codons. (in clone ". $fe->clone_id .")\n\n";
			next;
		    }
		    $tseq->desc("Gene:$gene_id.$gene_version Clone:".$fe->clone_id . " Contig:" . $fe->contig_id . " Chr: " . $chr . " Pos: " . $gene_start);
		    $seqiop->write_seq($tseq);
		    $seqiop=undef;
		    close(PEP);
		}
	    }
	    else {
		die "No valid format!";
	    }
	}
    };
    
    if( $@ ) {
	my $gene_id = "@chunk_list";
	#foreach my $clone ( $db->geneid_to_cloneid($gene_id)) {
	#    print ERR "Error in clone $clone:\n";
	#}
	print ERR "unable to process @chunk_list, due to \n$@\n";
    }
    $db->DESTROY;
    $db=undef;
}

sub  find_trans_start {
 my ($trans) = @_;

 my $start_pos;
 my $trans_pos;

 my $contig; 
 foreach my $exon ($trans->each_Exon) {
    if ($exon->id eq $trans->translation->start_exon_id) {
       $contig = $exon->contig_id;
       if ($exon->strand == 1) {
         $start_pos = $exon->start;
         $trans_pos = $exon->start + $trans->translation->start - 1;
       } else {
         $start_pos = $exon->end;
         $trans_pos = $exon->end - $trans->translation->start + 1;
       }
    }
 }
 if (!defined($start_pos)) {
    print ERR "Couldn't find start exon for " . $trans->id . "\n";
    die;
 }
 my $query = "select chr_name,chr_start,chr_end,raw_start,raw_end,raw_ori from static_golden_path,contig where raw_id = internal_id and contig.id = \'" . $contig . "\'";
 my $sth  = $db->prepare($query);
 my $res  = $sth->execute;

 my $row = $sth->fetchrow_hashref;

 my $chr = $row->{chr_name};
 my $chr_start = $row->{chr_start};
 my $chr_end   = $row->{chr_end};
 my $raw_start = $row->{raw_start};
 my $raw_end   = $row->{raw_end}; 
 my $raw_ori   = $row->{raw_ori};

 my $gene_start; 
 my $cdna_start;

 if ($raw_ori == 1) {
    $gene_start = $chr_start + ($trans_pos - $raw_start);
    $cdna_start = $chr_start + ($start_pos - $raw_start);
 } else {
    $cdna_start = $chr_end   - ($start_pos - $raw_start);
    $gene_start = $chr_end   - ($trans_pos - $raw_start);
 } 

 return ($chr,$gene_start,$cdna_start);
}


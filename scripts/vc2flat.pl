#!/usr/local/bin/perl 

=head1 NAME - vc2flat.pl

 provides flat file formats from Virtual Contigs in EnsEMBL databases

=head1 SYNOPSIS - 

    vc2flat chr1:10000-110000
  
    vc2flat -format fasta chr1:10000-11000

    vc2flat chr1

    vc2flat ctg1234

    vc2flat -format gff chr1

=head1 OPTIONS

    -verbose   print to STDERR on each clone to dump

    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::Obj)
    
    -host      host name for database (gets put as host= in locator)
			      
    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -nodna     dont write dna part of embl file (for testing)

    -format    [gff/ace/pep/fasta] dump in gff/ace/peptide/fasta format, not EMBL

    -pepformat What format output to dump to. NB gene2flat is a better
               script now for dumping translations.

    -usefile   read in on stdin a list of vc parameters, one per line

    -outfile   write output into file instead of to STDOUT

    -help      displays this documentation with PERLDOC

    -genetype  type of gene to get out of the database

=cut

use strict;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::EMBL_Dump;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::GTF_handler;

use Getopt::Long;

# global defaults
my $module    = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $dbtype    = 'rdb';
my $format    = 'embl';
my $nodna     = 0;
my $help;
my $noacc     = 0;
my $aceseq;
my $usefile  = 0;
my $pepformat = 'Fasta';
my $test=0;
my $part=0;
my $verbose=0;
my $outfile;
my $checkdna;
# test
my $species='';
my $freeze=0;
my $nogene=0;
my $nosecure=0;
my $genetype = undef;

# defaults for msql (rdb) access
my $host     = 'ensrv3';
my $dbname    = 'ensembl100';
my $dbuser    = 'ensro';
my $dbpass = undef;

# defaults for acedb (humace)
my $host2     = 'humsrv1';
my $port      = '410000';

my $focuscontig;
my $focusposition;
my $ori;
my $left;
my $right;

my $static = 0;

&GetOptions( 'module:s'  => \$module,
	     'dbtype:s'  => \$dbtype,
	     'dbuser:s'  => \$dbuser,
	     'dbpass:s'  => \$dbpass,
	     'host:s'    => \$host,
	     'port:n'    => \$port,
	     'dbname:s'  => \$dbname,
	     'format:s'  => \$format,
	     'nodna'     => \$nodna,
	     'h|help'    => \$help,
	     'noacc'     => \$noacc,
	     'aceseq:s'  => \$aceseq,
	     'pepform:s' => \$pepformat,
	     # usefile to be consistent with other scripts
	     'usefile'   => \$usefile,
	     'test'      => \$test,
	     'part'      => \$part,
	     'checkdna:s'=> \$checkdna,
	     'verbose'   => \$verbose,
	     'outfile:s' => \$outfile,
	     'species:s' => \$species,
	     'freeze:n'  => \$freeze,
	     'static'    => \$static,
	     'nogene'    => \$nogene,
	     'nosecure'  => \$nosecure,
	     'fcontig:s' => \$focuscontig,
	     'fposition:i' => \$focusposition,
	     'ori:i'     => \$ori,
	     'left:i'    => \$left,
	     'right:i'   => \$right,
	     'genetype:s' => \$genetype
	     ) or exec('perldoc', $0);

if ($help){
    exec('perldoc', $0);
}


# build database
my $db;
my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
$db = Bio::EnsEMBL::DBLoader->new($locator);
$db->static_golden_path_type('UCSC');
my $stadp = $db->get_StaticGoldenPathAdaptor();

#get input strings
my @vcstrings;
if( $usefile == 1 ) {
    while( <> ) {
	my ($string) = split;
	push(@vcstrings);
    }
} else {
    @vcstrings = @ARGV;
}

# set output file
my $OUT;
if($outfile){
    open(OUT,">$outfile") || die "cannot write to $outfile";
    $OUT=\*OUT;
}else{
    $OUT=\*STDOUT;
}


# sequence output files opened only once. Notice monkeying around
# with embl/genbank. Needs smarts to handle issues in there.
my $seqout;
if( $format =~ /embl|genbank/ ) {
    $seqout = Bio::SeqIO->new( '-format' => $format, -fh => $OUT);
    &Bio::EnsEMBL::EMBL_Dump::ensembl_annseq_output($seqout);
    if( $nodna == 1 ) {
	$seqout->_show_dna(0);
    }
} elsif ( $format =~ /fasta/ ) {
    $seqout = Bio::SeqIO->new( '-format' => 'Fasta' , -fh => $OUT);
} 

# Main loop over strings
foreach my $vcstring ( @vcstrings ) {

    eval {
	my $vc;
	if( $vcstring =~ /(chr\S+):(\d+)-(\d+)/ ) {
	    my $chr   = $1;
	    my $start = $2;
	    my $end   = $3;
	    $vc = $stadp->fetch_VirtualContig_by_chr_start_end($chr,$start,$end);
	} elsif ( $vcstring =~ /^(chr\S+)$/ ) {
	    my $chr   = $1;
	    $vc = $stadp->fetch_VirtualContig_by_chr_name($chr);
	} elsif ( $vcstring =~ /^(ctg\S+)$/ ) {
	    my $ctg   = $1;
	    $vc = $stadp->fetch_VirtualContig_by_fpc_name($ctg);
	}
	$vc->id($vcstring);
	if( $format =~ /gff/ ) {
	    my @seqfeatures = $vc->top_SeqFeatures();
	    foreach my $sf (@seqfeatures ) {
		print $OUT $sf->gff_string, "\n";
	    }
	} elsif ( $format =~ /fasta/ ) {
	    $seqout->write_seq($vc->primary_seq());
	} elsif ( $format =~ /embl/ ) {
	    &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($vc);
	    $seqout->write_seq($vc);
	} elsif ( $format =~ /genbank/ ) {
	    &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($vc);
	    $seqout->write_seq($vc);
	} elsif ( $format =~ /pep/ || $format =~ /gtf/ ) {
	    my @genes;
	    if( defined $genetype ) {
		@genes = $vc->get_Genes_by_Type($genetype);
	    } else {
		@genes = $vc->get_all_Genes();
	    }

	    my $vcid = $vc->id();

	    if( $format =~ /pep/ ) {
		foreach my $gene ( @genes ) {
		    my $geneid = $gene->id();
		    foreach my $trans ( $gene->each_Transcript ) {
			my $tseq = $trans->translate();
			if( $tseq->seq =~ /\*/ ) {	
			    print STDERR $trans-id," got stop codons. Not dumping. on $vcid\n";
			}
			$tseq->desc("VirtualContig:$vcid Gene:$geneid");
			$seqout->write_seq($tseq);
		    }
		}
	    } elsif ( $format =~ /gtf/ ) {
	      Bio::EnsEMBL::Utils::GTF_handler->dump_genes($OUT,@genes);
	    }
	} elsif ( $format =~ /ace/ ) {
	    $vc->write_acedb($OUT,$aceseq);
	}
    };
    if( $@ ) {
	print STDERR "Dumping Error! Cannot dump ".$vcstring."\n$@\n";
    }
}




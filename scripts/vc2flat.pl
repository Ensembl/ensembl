#!/usr/local/bin/perl 

=head1 NAME - vc2flat.pl

 provides flat file formats from Virtual Contigs in EnsEMBL databases

=head1 SYNOPSIS - 

    vc2flat -vc_parameters

    vc2flat -format gff 
   
    vc2flat -dbtype ace 

    vc2flat -dbtype rdb -host mysql.server.somewhere 

=head1 OPTIONS

    -verbose   print to STDERR on each clone to dump

    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::Obj)
    
    -species   Species other than human (e.g. 'mouse' - if not set, only loads human)

    -freeze    Loads only clones in current frozen set (number)

    -nogene    Allows clones to be read if dna has been updated successfully
               and searches completed, but before genes have been built and clone unlocked

    -nosecure

    -dbtype    Database type (valid types are rdb, ace)

    -host      host name for database (gets put as host= in locator)
			      
    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -nodna     dont write dna part of embl file (for testing)

    -format    [gff/ace/pep/fasta] dump in gff/ace/peptide/fasta format, not EMBL

    -pepformat What format output to dump to. NB gene2flat is a better
               script now for dumping translations.

    -usefile   read in on stdin a list of vc parameters, one set per line

    -outfile   write output into file instead of to STDOUT

    -help      displays this documentation with PERLDOC

=cut

use strict;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::EMBL_Dump;
use Bio::EnsEMBL::DB::VirtualContig;
use Bio::EnsEMBL::DB::EmblVirtualContig;
use Bio::SeqIO;

use Getopt::Long;

# signal handler
$SIG{INT}=sub {my $sig=shift;die "exited after SIG$sig";};

# global defaults
my $host = 'localhost';
my $module    = 'Bio::EnsEMBL::DBSQL::Obj';
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

# defaults for msql (rdb) access
my $host1     = 'localhost';
my $dbname    = 'ensembl08';
my $dbuser    = 'root';
my $dbpass = undef;

# defaults for acedb (humace)
my $host2     = 'humsrv1';
my $port      = '410000';

my $focuscontig;
my $focusposition;
my $ori;
my $left;
my $right;

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
	     'nogene'    => \$nogene,
	     'nosecure'  => \$nosecure,
	     'fcontig:s' => \$focuscontig,
	     'fposition:i' => \$focusposition,
	     'ori:i'     => \$ori,
	     'left:i'    => \$left,
	     'right:i'   => \$right
	     ) or exec('perldoc', $0);

if ($help){
    exec('perldoc', $0);
}

# set module if dbtype set and recognised


my $db;

my @vcs;

if( $usefile == 1 ) {
    while( <> ) {
	my ($focus_contig,$focus_position,$ori,$left,$right);
	if (/(.+)\,(.+)\,(.+)\,(.+)\,(.+)/) { 
	    $focus_contig=$1;
	    $focus_position=$2;
	    $ori=$3;
	    $left=$4;
	    $right=$5;
	}
	else {
	    print STDERR "Got $focus_contig and $focus_position and $ori and $left and $right\n";
	    die ("Could not read vc parameter file!\n The format should be one vc 
per line, ith vc parameters separated by commas!");
	}
	my @list=($focus_contig,$focus_position,$ori,$left,$right);
	push(@vcs,\@list);
    }
} else {
    $focuscontig || die ("Need to supply a focus contig!\n");
    $focusposition || die ("Need to supply a focus position!\n");
    $ori || die ("Need to supply an orientation!\n");
    $left || die ("Need to supply left size!\n"); 
    $right || die ("Need to supply right size!");
    
    my @list=($focuscontig,$focusposition,$ori,$left,$right);
    push(@vcs,\@list);
}


my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
$db = Bio::EnsEMBL::DBLoader->new($locator);

# set output file
my $OUT;
if($outfile){
    open(OUT,">$outfile") || die "cannot write to $outfile";
    $OUT=\*OUT;
}else{
    $OUT=\*STDOUT;
}


# Main loop over clone ids...
foreach my $vc_list_ref ( @vcs ) {

    if ($verbose) {
	print STDERR "Creating VirtualContig:\n";
	foreach my $parameter ($vc_list_ref) {
	    print STDERR "$parameter ";
	}
	print STDERR "\n";
    }

    eval {
	my $focuscontig=$db->get_Contig($vc_list_ref->[0]);
	my $vc=Bio::EnsEMBL::DB::EmblVirtualContig->new(
	 				       -focuscontig => $focuscontig,
	       				       -focusposition => $vc_list_ref->[1],
	       				       -ori => $vc_list_ref->[2],
	       				       -left => $vc_list_ref->[3], 
	       				       -right => $vc_list_ref->[4] 
						       );
	$vc->sv(1);
        # debug tests by contig
	print(STDERR "Format is $format\n");
	if( $format =~ /gff/ ) {
	    my @seqfeatures = $vc->top_SeqFeatures();
	    foreach my $sf (@seqfeatures ) {
		print $OUT $sf->gff_string, "\n";
	    }
	} elsif ( $format =~ /fasta/ ) {
	    my $seqout = Bio::SeqIO->new( '-format' => 'Fasta' , -fh => $OUT);
	    $seqout->write_seq($vc->primary_seq());
	} elsif ( $format =~ /embl/ ) {
	    print(STDERR "Dumping embl\n");
	    &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($vc);
	    print(STDERR "making new seq\n");
	    my $emblout = Bio::SeqIO->new( '-format' => 'EMBL', -fh => $OUT);
	    &Bio::EnsEMBL::EMBL_Dump::ensembl_annseq_output($emblout);
	    if( $nodna == 1 ) {
		$emblout->_show_dna(0);
	    }
	    $emblout->write_seq($vc);
	} elsif ( $format =~ /genbank/ ) {
	    &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($vc);
	    my $gbout = Bio::SeqIO->new( '-format' => 'GenBank', -fh => $OUT);
	    &Bio::EnsEMBL::EMBL_Dump::ensembl_annseq_output($gbout);
	    # genbank format - the ID line is wrong. Fall back to locus
	    $gbout->_id_generation_func(undef);
	    $gbout->_ac_generation_func(undef);
	    $vc->accession($vc->id());
	    $vc->division("PRI");
	    if( $nodna == 1 ) {
		$gbout->_show_dna(0);
	    }
	    $gbout->write_seq($vc);
	} elsif ( $format =~ /pep/ ) {
	    my $seqout = Bio::SeqIO->new ( '-format' => $pepformat , -fh => $OUT ) ;
	    my $vcid = $vc->id();
	    foreach my $gene ( $vc->get_all_Genes() ) {
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
	} elsif ( $format =~ /ace/ ) {
	    $vc->write_acedb($OUT,$aceseq);
	}
    };
    if( $@ ) {
	print STDERR "Dumping Error! Cannot dump ".$vc_list_ref->[0]."\n$@\n";
    }
}



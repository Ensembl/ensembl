#!/usr/local/bin/perl -w

=head1 NAME - cam2flat.pl

 provides flat file format from Camace databases

=head1 SYNOPSIS - 

    cam2flat dJ271M21

    cam2flat -format gff dJ271M21
   
    cam2flat -dbtype ace dJ271M21

    cam2flat -dbtype rdb -host mysql.server.somewhere dJ271M21


=head1 OPTIONS

    -verbose   print to STDERR on each clone to dump

    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::Obj)
    
    -dbtype    Database type (valid types are timdb, ace)

    -host    host name for database (gets put as host= in locator)
			      
    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -nodna     dont write dna part of embl file (for testing)

    -format    [gff/ace/pep/fasta] dump in gff/ace/peptide/fasta format, not EMBL

    -pepformat What format output to dump to. NB gene2flat is a better
               script now for dumping translations.

    -noacc     [only timdb] by default, regardless of specifing the
               accession for a sanger clone or its clonename, it will
               dump as its accession.  Use -noacc to dump by clonename

    -getall    all clones from the database [no applicable to timdb]

    -start     start point in list of clones (useful with -getall)

    -end       end point in list of clones (useful with -getall)

    -outfile   write output into file instead of to STDOUT

    -help      displays this documentation with PERLDOC


=cut



BEGIN
{
  my $rootdir = "/nfs/disk65/sjk/perl";
  unshift (@INC,"$rootdir/ensembl/modules");
  unshift (@INC,"$rootdir/PerlModules");
}

use lib qw(/nfs/disk100/humpub/modules/bioperl-0.6);

use strict;
use Bio::EnsEMBL::AceDB::Contig;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::EMBL_Dump;
use Bio::SeqIO;
use Getopt::Long;

# signal handler
$SIG{INT}=sub {my $sig=shift;die "exited after SIG$sig";};

# global defaults
my $host      = 'wormsrv1';
my $module    = 'Bio::EnsEMBL::AceDB::Obj';
my $dbtype    = '';
my $format    = 'embl';
my $nodna     = 0;
my $help;
my $noacc     = 0;
my $aceseq;
my $getall    = 0;
my $pepformat = 'Fasta';
my $part;
my $live;
my $verbose   = 0;
my $cstart    = 0;
my $cend      = undef;
my $outfile;
my $oldstyle = 0;
my $checkdna;
my $host1     = 'obi-wan';
my $dbname    = 'ensdev';
my $dbuser    = 'ensembl';
my $dbpass = 'pass';
# defaults for acedb (humace)
my $host2     = 'humsrv1';
my $port      = '100100';


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
	     'getall'    => \$getall,
	     'part'      => \$part,
	     'live'      => \$live,
	     'checkdna:s'=> \$checkdna,
	     'verbose'   => \$verbose,
	     'outfile:s' => \$outfile,
	     'oldstyle'  => \$oldstyle,
	     ) or exec('perldoc', $0);




if ($help) {
    exec('perldoc', $0);
}

my @clones = @ARGV;
my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db = Bio::EnsEMBL::DBLoader->new($locator);

# Note because of a bug with this version of AcePerl two clones F47A4 and W05E10 do not get retreived
if ( $getall == 1 ) {
    @clones = $db->get_all_Clone_id();    
    print STDERR scalar(@clones)." clones found in DB\n";
}

# set output file
my $OUT;
if($outfile){
    open(OUT,">$outfile") || die "cannot write to $outfile";
    $OUT=\*OUT;
}else{
    $OUT=\*STDOUT;
}


#
# Main loop over clone ids...
#
#
#

foreach my $clone_id ( @clones ) {

    if( $verbose >= 1 ) {
	print STDERR "Dumping $clone_id\n";
    }

    #
    # wrap in eval to protect from exceptions being
    # thrown. One problem here is that if the exception
    # is thrown during the feature table dumping, then it
    # has still dumped part of the entry. Bad news. Need to "exercise" the
    # objects first.
    #
    # I think we just have to assumme that if we get an object that
    # it conforms to the interface. 
    #

    eval {
	my $clone = $db->get_Clone($clone_id);
	my $as    = $clone->virtualcontig;
	$as->skip_SeqFeature('similarity',1);
        
	# choose output mode
	
	# test clone check call
	if($checkdna){
	    if($clone->compare_dna($checkdna)){
		print STDERR "WARN: DNA of $checkdna DIFFERENT to that in $clone_id\n";
	    }else{
		print STDERR "DNA of $checkdna identical to that in $clone_id\n";
	    }
	}

	my $debug = 0;
	if($debug){
	    # debug tests
	    my $id=$clone->id;
	    print "DNA date [$id]: ".$clone->seq_date."\n";
	    print "create date [$id]: ".$clone->created."\n";
	    print "modify date [$id]: ".$clone->modified."\n";
	    print "version [$id]: ".$clone->version."\n";

	    # debug tests by contig
	    foreach my $contig ($clone->get_all_Contigs){
		my $id=$contig->id;
		print "DNA date [$id]: ".$contig->seq_date."\n";
	    }
	    exit 0;
	}
        
        
	if( $format =~ /gff/ ) {
	    foreach my $contig ( $clone->get_all_Contigs )  {
		my @seqfeatures = $contig->as_seqfeatures();
		foreach my $sf ( @seqfeatures ) {
		    print $OUT $sf->gff_string, "\n";
		}
	    }
	} 
        
        elsif ( $format =~ /test/ ) {
	    foreach my $contig ( $clone->get_all_Contigs() ) {
            
                my  $vc = Bio::EnsEMBL::DB::VirtualContig->new( 
                                                -focuscontig => $contig,
                                                  -focusposition => 10000,
                                                  -ori => 1,
                                                  -left => 80000,
                                                  -right => 70000
                                                  );            
                print STDERR "Created Contig, lenght = ", $vc->length(), "\n";
                my $seq = $vc->primary_seq();
                print STDERR "primary sequence lenght = ", $seq->length(), "\n";
            }
        }
            
            
        elsif ( $format =~ /fasta/ ) {
	    my $seqout = Bio::SeqIO->new( '-format' => 'Fasta' , -fh => $OUT);
	    
	    foreach my $contig ( $clone->get_all_Contigs() ) {
                
                # Create a virtual contig around the contig
                my  $vc = Bio::EnsEMBL::DB::VirtualContig->new( 
                                                -focuscontig => $contig,
                                                  -focusposition => 1,
                                                  -ori => 1,
                                                  -left => 100000,
                                                  -right => 100000
                                                  );            
                print STDERR "Created Virtual Contig, lenght = ", $vc->length(), "\n";
                my $seq = $vc->primary_seq();
                print STDERR "primary sequence lenght = ", $seq->length(), "\n";
                
                # Loop through all the genes on the virtual contig
                for my $gene ($vc->get_all_Genes()) {
          
                    # Loop through all the transcripts of each gene
                    for my $transcript ($gene->each_Transcript()) {

                        # Get the translation 
                        my $translation = $transcript->translate();                        
                        # Check that the translation does not contain '*'s                        
                        if ($translation =~ m/\*/) {
                            print STDERR "* found in sequence of $gene->id()";
                        }  
                        print $OUT "Gene ID: ", $gene->id, "\n";
                        print $OUT "First Exon: ", $transcript->first_exon()->id, "\n";
                        print $OUT "Last Exon: ", $transcript->last_exon()->id, "\n";
                        print $OUT "Start Exon: ", $transcript->start_exon()->id, "\n";
                        print $OUT "End Exon: ", $transcript->end_exon()->id, "\n";                                                 
                        print $OUT "Translation: ", $translation->seq(), "\n";
                        print STDERR "Written translation of ", $gene->id, "\n";
                    }
                }
	    }
	} 
        
    
        elsif ( $format =~ /embl/ ) {
	    print(STDERR "Dumping embl\n");
	    &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($as);
	    print(STDERR "making new seq\n");
	    my $emblout = Bio::SeqIO->new( '-format' => 'EMBL', -fh => $OUT);
	    &Bio::EnsEMBL::EMBL_Dump::ensembl_annseq_output($emblout);
	    if( $nodna == 1 ) {
		$emblout->_show_dna(0);
	    }
	    
	    $emblout->write_seq($as);
	} 
        
        elsif ( $format =~ /genbank/ ) {
	    &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($as);
	    my $gbout = Bio::SeqIO->new( '-format' => 'GenBank', -fh => $OUT);
	    &Bio::EnsEMBL::EMBL_Dump::ensembl_annseq_output($gbout);
	    # genbank format - the ID line is wrong. Fall back to locus
	    $gbout->_id_generation_func(undef);
	    $gbout->_ac_generation_func(undef);
	    $as->accession($clone->id());
	    $as->division("PRI");
	    if( $nodna == 1 ) {
		$gbout->_show_dna(0);
	    }
	    
	    $gbout->write_seq($as);
	} 
        
        elsif ( $format =~ /pep/ ) {
	    my $seqout = Bio::SeqIO->new ( '-format' => $pepformat , -fh => $OUT ) ;
	    my $cid = $clone->id();
	    
	    foreach my $gene ( $clone->get_all_Genes() ) {
		my $geneid = $gene->id();
		foreach my $trans ( $gene->each_Transcript ) {
		    my $tseq = $trans->translate();
		    if( $tseq->seq =~ /\*/ ) {	
			print STDERR $trans-id," got stop codons. Not dumping. on $cid\n";
		    }
		    $tseq->desc("Clone:$cid Gene:$geneid");
		    $seqout->write_seq($tseq);
		}
	    }
	} 
        
        elsif ( $format =~ /ace/ ) {
	    foreach my $contig ( $clone->get_all_Contigs() ) {
		$contig->write_acedb($OUT,$aceseq);
	    }
	}
    };
    
    if( $@ ) {
	print STDERR "Dumping Error! Cannot dump $clone_id\n$@\n";
    }
}



#!/usr/local/bin/perl

=head1 NAME - clone2embl

 provides flat file formats from EnsEMBL databases

=head1 SYNOPSIS - 

    clone2embl dJ271M21

    clone2embl -gff dJ271M21
   
    clone2embl -dbtype ace dJ271M21

    clone2embl -dbtype rdb -host mysql.server.somewhere dJ271M21

    clone2embl -dbtype timdb AL035541           # dump as accession
    clone2embl -dbtype timdb dJ718J7            # dump as accession
    clone2embl -dbtype timdb -noacc dJ718J7     # dump as clone

=head1 OPTIONS

    -dbtype    database type (rdb, ace, timdb)

    -nodna     don't write dna part of embl file (for testing)

    -format    [gff/ace/pep] dump in gff/ace/peptides format, not EMBL

    -noacc     by default, regardless of specifing the accession for a sanger clone 
               or its clonename, it will dump as its accession.  Use -noacc to 
               dump by clonename

    -test      use test database rather than live [currently only timdb]
               clones in testdb are listed with a T below

=head1 EXAMPLE CLONES

    dJ271M21/AL031983   T  single contig, mainly forward strand genes, but one reverse

    dJ718J7                single contig, single reverse strand gene with partial transcripts

    C361A3/AL031717     T  unfinished sanger, 3 contigs, 2 ordered by gene predictions

    C367G8/Z97634       T  unfinished sanger, 2 contigs, not ordered by gene predictions

    AP000228               External finished clone

    AC005776               External finished clone

    AC010144               External unfinished clone

=cut

# to find more things in TimDB use:
#~humpub/scripts.devel/set_dbm.pl -f ~/th/unfinished_ana/unfinished_clone -l 10

use strict;

#use Bio::EnsEMBL::AceDB::Obj;
use Bio::EnsEMBL::DB::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::AnnSeqIO;
use Bio::SeqIO;

use Getopt::Long;

my $dbtype = 'rdb';
my $host;
my $host1  = 'croc';
my $host2  = 'humsrv1';
my $port   = '410000';
my $format = 'embl';
my $nodna  = 0;
my $help;
my $noacc  = 0;
my $aceseq;
my $fromfile = 0;
my $getall =0;
my $pepformat = 'Fasta';
my $test;

# this doesn't have genes (finished)
#my $clone  = 'dJ1156N12';
# this does have genes (finished)
#my $clone  = 'dJ271M21';
# this does have genes (unfinished)
# my $clone = '217N14';

&GetOptions( 'dbtype:s'  => \$dbtype,
	     'host:s'    => \$host,
	     'port:n'    => \$port,
	     'format:s'  => \$format,
	     'nodna'     => \$nodna,
	     'h|help'    => \$help,
	     'noacc'     => \$noacc,
	     'aceseq:s'  => \$aceseq,
	     'pepform:s' => \$pepformat,
	     'fromfile'  => \$fromfile,
	     'getall'    => \$getall,
	     'test'      => \$test,
	     );

if($help){
    exec('perldoc', $0);
}

my $db;



if( $dbtype =~ 'ace' ) {
    $host=$host2 unless $host;
    $db = Bio::EnsEMBL::AceDB::Obj->new( -host => $host, -port => $port);
} elsif ( $dbtype =~ 'rdb' ) {
    $host=$host1 unless $host;
    $db = Bio::EnsEMBL::DB::Obj->new( -user => 'root', -db => 'ensdev' , -host => $host );
} elsif ( $dbtype =~ 'timdb' ) {

    # EWAN: no more - you should be able to load as many clones as you like!
    if( $#ARGV == -1 ) {
	push(@ARGV,'dJ271M21');
    }

    # clones required are passed to speed things up - cuts down on parsing of flat files
    $db = Bio::EnsEMBL::TimDB::Obj->new(\@ARGV,$noacc,$test);
} else {
    die("$dbtype is not a good type (should be ace, rdb or timdb)");
}

my @clones;

if( $fromfile == 1 ) {
    while( <> ) {
	my ($cloneid) = split;
	push(@clones,$cloneid);
    }
} elsif ( $getall == 1 ) {
    @clones = $db->get_all_Clone_id();
} else {
    @clones = @ARGV;
}


foreach my $clone_id ( @clones ) {

    eval {
	my $clone = $db->get_Clone($clone_id);
	my $as = $clone->get_AnnSeq();
	# choose output mode
	
	
	if( $format =~ /gff/ ) {
	    foreach my $contig ( $clone->get_all_Contigs )  {
		my @seqfeatures = $contig->as_seqfeatures();
		foreach my $sf ( @seqfeatures ) {
		    print $sf->gff_string, "\n";
		}
	    }
	} elsif ( $format =~ /fastac/ ) {
	    my $seqout = Bio::SeqIO->new( -format => 'Fasta' , -fh => \*STDOUT);
	    
	    foreach my $contig ( $clone->get_all_Contigs ) {
		$seqout->write_seq($contig->seq());
	    }
	} elsif ( $format =~ /embl/ ) {
	    
	    $as->seq->desc("Reannotated sequence via EnsEMBL");
	    my $comment = Bio::Annotation::Comment->new();
	    
	    $comment->text("This sequence was reannotated via the EnsEMBL system. Please visit the EnsEMBL web site, http://ensembl.ebi.ac.uk for more information");
	    $as->annotation->add_Comment($comment);
	    
	    $comment = Bio::Annotation::Comment->new();
	    $comment->text("The /gene_id indicates a unique id for a gene, /transcript_id a unique id for a transcript and a /exon_id a unique id for an exon. These ids are maintained wherever possible between versions. For more information on how to interpret the feature table, please visit http://ensembl.ebi.ac.uk/docs/embl.html");
	    $as->annotation->add_Comment($comment);
	    
	    my $sf = Bio::SeqFeature::Generic->new();
	    $sf->start(1);
	    $sf->end($as->seq->seq_len());
	    $sf->primary_tag('source');
	    $sf->add_tag_value('organism','Homo sapiens');
	    $as->add_SeqFeature($sf);
	    my $emblout = Bio::AnnSeqIO->new( -format => 'EMBL', -fh => \*STDOUT);
	    $emblout->_post_sort(\&sort_FTHelper_EnsEMBL);
	    
	    # attach ensembl specific dumping functions
	    $emblout->_id_generation_func(\&id_EnsEMBL);
	    $emblout->_kw_generation_func(\&kw_EnsEMBL);
	    $emblout->_sv_generation_func(\&sv_EnsEMBL);
	    $emblout->_ac_generation_func(\&ac_EnsEMBL);
	    
	    if( $nodna == 1 ) {
		$emblout->_show_dna(0);
	    }
	    
	    $emblout->write_annseq($as);
	} elsif ( $format =~ /pep/ ) {
	    my $seqout = Bio::SeqIO->new ( '-format' => $pepformat , -fh => \*STDOUT ) ;
	    
	    foreach my $gene ( $clone->get_all_Genes() ) {
		foreach my $trans ( $gene->each_Transcript ) {
		    my $tseq = $trans->translate();
		    $seqout->write_seq($tseq);
		}
	    }
	} elsif ( $format =~ /ace/ ) {
	    foreach my $contig ( $clone->get_all_Contigs() ) {
		$contig->write_acedb(\*STDOUT,$aceseq);
	    }
	}
    };
    if( $@ ) {
	print STDERR "Dumping Error! Cannot dump $clone_id\n$@\n";
    }
}
    
#########################
# sub routines
#########################

sub id_EnsEMBL {
    my $annseq = shift;

    return sprintf("%-11s standard; DNA; %s; %d BP.",$annseq->embl_id(),$annseq->htg_phase == 4 ? 'HUM' : 'HTG',$annseq->seq->seq_len() );
}


sub kw_EnsEMBL{
   my ($annseq) = @_;

   if( $annseq->htg_phase == 4 ) {
       return "HTG.";
   }

   return "HTG; HTG_PHASE" . $annseq->htg_phase() . ".";
}

sub sv_EnsEMBL {
   my ($annseq) = @_;

   if( ! $annseq->sv ) {
       return "NO_SV_NUMBER";
   }
   if( $annseq->sv == -1 ) {
       return undef;
   }

   return $annseq->seq->id() . "." . $annseq->sv
}

sub ac_EnsEMBL {
   my ($annseq) = @_;

   return $annseq->seq->id() . ";";
}


sub sort_FTHelper_EnsEMBL {
    my $a = shift;
    my $b = shift;

    if( $a->key eq $b->key ) {
	return ($a->loc cmp $b->loc);
    }

    if( $a->key eq 'source' ) {
	return -1;
    }
    if( $b->key eq 'source' ) {
	return 1;
    }

    if( $a->key eq 'CDS' ) {
	return -1;
    }

    return 1;
}


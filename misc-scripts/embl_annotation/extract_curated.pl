#!/usr/local/bin/perl

# script to extract annotation from EMBL files and insert it into an
# ensembl database

use strict;

use Getopt::Long;

use Bio::EnsEMBL::EMBLLOAD::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::SeqIO;

my $ohost='ecs1b';
my $ouser='ensadmin';
my $odbname='embl100';

my $ihost='ecs1c';
my $iuser='ensadmin';
my $idbname='ensembl100';

my $emblfile='AL109928.embl';
my $write;
my $verbose;
my $list;
my $delete;
my $clone_id;
my $processed_clone_list='extract_curated.processed.lis';
my $log='extract_curated.log';
my $max_clone;

&GetOptions(

	    'ohost:s'=>\$ohost,
	    'ouser:s'=>\$ouser,
	    'odbname:s'=>\$odbname,

	    'ihost:s'=>\$ihost,
	    'iuser:s'=>\$iuser,
	    'idbname:s'=>\$idbname,

	    'emblfile:s'=>\$emblfile,

	    'w'=>\$write,
	    'v'=>\$verbose,
	    'l'=>\$list,
	    'd'=>\$delete,
	    'c:s'=>\$clone_id,
	    'log:s'=>\$log,
	    'processed:s'=>\$processed_clone_list,
	    'max:n'=>\$max_clone,
	     );

# use -d -c clone to force a delete

my $nunfin=0;
my $nnoembl=0;
my $nwrongsv=0;
my $nclonefail=0;
my $nclonenogene=0;
my $nclonegene=0;
my $ngeneok=0;

my %processed;
my $logging;

# connect to db for writing
my $db = Bio::EnsEMBL::DBLoader->new(
    "Bio::EnsEMBL::DBSQL::DBAdaptor/host=$ohost;user=$ouser;dbname=$odbname");

my $dbi;
if($list){

    # open log
    if($write){
	open(LOG,">$log") || die "cannot open $log";
	$logging=1;
    }

    # load current list of processed clones
    if(open(IN,$processed_clone_list)){
	my $n=0;
	while(<IN>){
	    if(/^(\S+)\s+(\S+)/){
		$processed{$1}=$2;
		$n++;
	    }
	}
	close(IN);
	print "$n processed clones read from $processed_clone_list\n";
    }else{
	print "WARN: no $processed_clone_list file found\n";
    }

    $dbi=Bio::EnsEMBL::DBLoader->new(
    "Bio::EnsEMBL::DBSQL::DBAdaptor/host=$ihost;user=$iuser;dbname=$idbname");

    my $nclone=0;
    foreach my $clone_id ($dbi->get_all_Clone_id){

	# has clone already been processed?
	if($processed{$clone_id}){
	    next;
	}

	$nclone++;
	last if($max_clone && $nclone>$max_clone);

	# is clone in one contig (finished)?
	my @contig;
	my $clone=$dbi->get_Clone($clone_id);
	eval{
	    (@contig)=$clone->get_all_Contigs;
	};
	if($@){
	    print "Contig error in Database for $clone_id:\n$@\n\n";
	    next;
	}elsif(scalar(@contig)>1){
	    print "$clone_id is in ".scalar(@contig)." contigs - skip\n";
	    $processed{$clone_id}='unfinished';
	    $nunfin++;
	    next;
	}

	my $version=$clone->embl_version;

	# is clone in output database?
	eval{
	    my $clone=$db->get_Clone($clone_id);
	};
	if($@){
	    print "[$nclone] $clone_id not yet parsed ".
		"[$nunfin; $nnoembl; $nwrongsv; $nclonefail; ".
		    "$nclonenogene; $nclonegene; $ngeneok]\n";
	    my $seqio;
	    eval{
		$seqio = Bio::SeqIO->new(-format => 'EMBL', -file => 
		    "/usr/local/pubseq/bin/efetch emnew:$clone_id |");
	    };
	    if($@){
		print "$clone_id could not be fetched from EMBL - warn\n";
		$processed{$clone_id}='embl_failed';
		$nnoembl++;
	    }
	    &_process_file($db,$seqio,$verbose,$clone_id,$version);
	}else{
	    print "$clone_id already parsed - skip\n";
	}

    }

    # write current list of checked, empty clones
    if($write){
	rename($processed_clone_list,"$processed_clone_list.old");
	if(open(OUT,">$processed_clone_list")){
	    foreach my $clone_id (keys %processed){
		print OUT "$clone_id ".$processed{$clone_id}."\n";
	    }
	    close(OUT);
	}else{
	    print "ERROR: failed writing to $processed_clone_list file\n";
	}
	close(LOG);
    }

}else{

    if(!$clone_id && $emblfile=~/^(\w+)/){
	$clone_id=$1;
    }

    # delete this clone
    if($delete){

	# remove all genes
	foreach my $gene_id ($db->gene_Obj->get_all_Gene_id()){
	    $db->gene_Obj->delete($gene_id);
	}
	
	my $clone=$db->get_Clone($clone_id);

	# remove clone
	$clone->delete_by_dbID;
	exit 0;
    }

    # link to file to be read
    my $seqio = Bio::SeqIO->new(-format => 'EMBL',-file => $emblfile);
    &_process_file($db,$seqio,$verbose,$clone_id,0);

}

# loop over entries in file
sub _process_file{
    my($db,$seqio,$verbose,$clone_id,$version)=@_;

    for(;;){

	my $seq;
	eval{
	    $seq = $seqio->next_seq;
	};
	if($@){
	    print "Problem accessing EMBL file for $clone_id\n$@\n\n";
	    next;
	}
	if(!$seq){
	    last;
	}

	# check SV
	my $sv=$seq->sv;
	if($sv=~/\.(\d+)$/){
	    if($version && $1!=$version){
		print "Sequence versions are different: Ensembl:$version; EMBL:$1\n";
		$processed{$clone_id}="embl_version:ENSEMBL:$version:EMBL:$1";
		$nwrongsv++;
		next;
	    }else{
		print "Sequence versions are same: $version\n";
	    }
	}else{
	    print "SV could not be parsed |$sv|\n";
	    next;
	}

	# get sequence from embl entry
	my $obj = Bio::EnsEMBL::EMBLLOAD::Obj->new(-seq => $seq);

	# HACK HACK HACK
	# there is only one clone to get, so don't need an id, but
	# need to set the id afterwards since it is wrong (embl_id not accession)
	my $clone = $obj->get_Clone();
	$clone->id($clone_id);
	
	# only for finished sequence, so only one contig
	my $contig;
	($contig) = $clone->get_all_Contigs();

	my @genes;
	eval{
	    @genes=($contig->get_all_Genes);
	};
	if($@){
	    print "Problems processing clone $clone_id\n";
	    if($logging){
		print LOG "$clone_id\n$@\n\n";
	    }else{
		print "$@\n\n";
	    }
	    $nclonefail++;
	    next;
	}

	# loop over entries to get genes and insert them into database
	my $ngene=scalar(@genes);
	if($ngene){
	    $nclonegene++;

	    # only write clone entry into database if there are genes
	    if($write){
	      eval{
		$db->write_Clone($clone);
	      };
	      if($@){
		print "Failed to write clone $clone_id\n";
	      }
	    }

	}else{
	    # processing ok, but no genes
	    $processed{$clone_id}='no_genes';
	    $nclonenogene++;
	}
	$ngeneok+=$ngene;
	print "parsed $ngene ok [$ngeneok]\n";

	foreach my $gene ( @genes ) {
	    
	    # writeout information about gene
	    if($verbose){
		my $gene_id=$gene->id;
		print "$gene_id\n";
		foreach my $transcript ($gene->each_Transcript){
		    my $transcript_id=$transcript->id;
		    print "  $transcript_id\n";
		}
	    }
	    
	    # write gene to database
	    if($write){
	      eval{
		$db->write_Gene($gene);
	      };
	      if($@){
		print "Failed to write gene ".$gene->id."\n";
	      }
	    }
	}
    }
}





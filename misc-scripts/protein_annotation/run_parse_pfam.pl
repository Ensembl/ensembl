#!/ebi/supported/bin/perl -w
use strict;
use Getopt::Long;
use DBI;
use Bio::Tools::BPlite;
use IO::Handle;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::Feature_Obj;
use Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor;

my ($query,$sp,$dbname,$host,$dir);

my $script;

&GetOptions(
            
            'query:s'=>\$query,
	    'script:s'=>\$script,
            'data_dir:s'=>\$dir,
	    'dbname:s'=>\$dbname,
	    'host:s'=>\$host
            );

print STDERR "DIR: $script\n";

#&run_pfam();
my ($db) = &db_connect();
my ($analysis) = &get_write_analysis($db);
&split_blast($db,$analysis);
#disconnection($db);

sub run_pfam {
    print STDERR "DIR: $script\n";
    
    my $pfam= "csh $script/run_pfam.csh $dir/$query";
    system($pfam) == 0 or die "$0\Error running '$dir' : $!";
   
}

sub db_connect {
    print STDERR "Connecting to the database...\n";
    my $enslocator = "Bio::EnsEMBL::DBSQL::Obj/host=ecs1b;dbname=prot_pipeline_test;user=root;perlonlyfeatures=1";
    
    my $db =  Bio::EnsEMBL::DBLoader->new($enslocator);
    
    return ($db);
}

sub get_write_analysis {
    print STDERR "Get write Analysis\n";    

    my ($db)=@_;

    my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
    
    my $anal2    = new Bio::EnsEMBL::Analysis(-db              => "Pfam",
					      -db_version      => 1,
					      -program         => "hmmpfam",
					      -program_version => 1,
					      -gff_source      => "hmmpfam",
					      -gff_feature     => "similarity",
					      );
#write the analysis into the database.
    $feature_obj->write_Analysis($anal2);
    

#my $query = $db->prepare("insert into analysis values(NULL,'Pfam',1,'hmmpfam',1,'hmmpfam','similarity');");
			
	    #$query->execute;
    

    return $anal2;
}


sub split_blast {
    print STDERR "Splitting and loading blast\n";
    my ($db,$analysis) = @_;
    my $count = 1;
    open (README1,"$query.out") || die "cant open README";
    
    while (<README1>) {
	if (($_ =~ /^HMM/) && ($count == 1)) {
	    
	    open (README, ">$dir/out.1");
	    $count++;
	}
	elsif($_ !~ /^HMM/) {
	    print README "$_";
	}
	elsif (($_ =~ /^HMM/) && ($count != 1)) {
	    close (README);
	    &parse_blast($db,$analysis);
	    my $rm = "rm $dir/out.1";
	    system($rm)==0 || die "$0\Error running '$rm' : $!";
	    open (README, ">out.1");
	    $count++;
	}
	
    }
}

sub parse_blast {
    my ($db,$analysis) = @_;
    
    open (README2,"$dir/out.1") || die "cant open README";
    
#Use BPlite to parse the Pfam output which is in blast format (more or less)   
    my $report = new Bio::Tools::BPlite(-fh=>\*README2);

    #Use protein feature adaptor to write the protein features in the database  
    my $protfeat = Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor->new($db);
        
    while(my $sbjct = $report->nextSbjct) {
	
	while (my $hsp = $sbjct->nextHSP) {
	    
	    
#Put the different components of the Pfam output to make it clear	    
	    my($ensemblAC) = $report->query =~ /\A(\w+)/;	
	    
	    my $idt = (int ($hsp->percent + 0.5));
	    my ($pfam) = $sbjct->name =~ /\A(\w+)/;
	    
	    my $seq_start = $hsp->query->start;
	    my $seq_end = $hsp->query->end;
	    my $subj_start = $hsp->subject->start;
	    my $subj_end = $hsp->subject->end;
	    my $transl = "transl";
	    my $anal = 22;
	    my $hid = $pfam;
	    my $score = int($hsp->score);
	    my $evalue = $hsp->P;
	    my $perc_id = $hsp->percent;
	    
	    if ($evalue <= 0.01) {
		
		my $feat1 = new Bio::EnsEMBL::SeqFeature ( -start => $seq_start,                   
							   -end => $seq_end,        
							   -score => $score,
							   -analysis => $analysis,
							   -seqname => "$ensemblAC",
							   -percent_id => $perc_id,
							   -p_value => "$evalue");
		
		my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => $subj_start,
							  -end => $subj_end,
							  -analysis => $analysis,
							  -seqname => "$pfam");
		
		
		my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
							    -feature2 => $feat2);
		
		$protfeat->write_Protein_feature($feature);
	    }
	}
    }
}



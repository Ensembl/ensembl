use strict;

use lib '../../modules/','../../../bioperl-live';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Lite::DBAdaptor;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname  );


GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname
	  );
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $user,
					    -port => $port,
					    -pass => $pass,
					    -dbname => $dbname);

#
# Only run on database with genes
#
my $sth = $db->prepare( "select count(*) from gene" );
$sth->execute();

my ( $gene_count )  = $sth->fetchrow_array();

if( ! $gene_count ) {
    print STDERR "No gene density for $dbname.\n";
    exit();
}

#
# and seq_regions
#
$sth = $db->prepare( "select count(*) from seq_region" );
$sth->execute();
my ( $seq_region_count ) = $sth->fetchrow_array();
if( ! $seq_region_count ) {
    print STDERR "No seq_regions for $dbname.\n";
    exit();
}


#my $snps_present = lite_attach( $db );
my $snps_present = 0;


my $slice_adaptor = $db->get_SliceAdaptor();
my $attrib_adaptor = $db->get_AttributeAdaptor();

my $top_slices = $slice_adaptor->fetch_all( "toplevel" );


foreach my $slice (@$top_slices) {
    my $num_known_genes  = 0;
    my $num_nov_CDS      = 0;
    my $num_nov_trans    = 0;
    my $num_tot_pseudo_genes = 0;
    my $num_unclass_pseudo_genes = 0;
    my $num_proc_pseudo_genes = 0;
    my $num_unproc_pseudo_genes = 0;
    my $num_put_trans = 0;
    my $num_pred_trans = 0;
    my $num_pred_Ig_pseudogenes = 0;
    my $num_Ig_segments = 0;
    
    print STDERR "Processing seq_region ", $slice->seq_region_name(), "\n";

    my @genes = @{$slice->get_all_Genes()};

    foreach my $gene (@genes) {
	my $type = $gene->type;
	if ($type eq 'Known') {
	    $num_known_genes++;
	} elsif ($type eq 'Novel_CDS') {
	    $num_nov_CDS++;
	} elsif ($type eq 'Novel_Transcript') {
	    $num_nov_trans++;
	} elsif ($type eq 'Putative') {
	    $num_put_trans++;
	} elsif ($type eq 'Predicted_Gene') {
	    $num_pred_trans++;
	} elsif ($type eq 'Ig_Pseudogene_Segment') {
	    $num_pred_Ig_pseudogenes++;
	} elsif ($type eq 'Ig_Segment') {
	    $num_Ig_segments++;
	} elsif ($type eq 'Pseudogene') {
	    $num_unclass_pseudo_genes++;
	} elsif ($type eq 'Processed_pseudogene') {
	    $num_proc_pseudo_genes++;
	} elsif ($type eq 'Unprocessed_pseudogene') {
	    $num_unproc_pseudo_genes++;
	}

	if ($type =~ /seudogene$/) {
	    $num_tot_pseudo_genes++;
	}
    }

    print "Slice", $slice->seq_region_name(), " has the following features:\n\n";
    print "known genes = $num_known_genes\n";
    print "novel coding sequences = $num_nov_CDS\n";
    print "novel transcripts = $num_nov_trans\n";
    print "putative transcripts = $num_put_trans\n";
    print "predicted transcripts = $num_pred_trans\n"; 
    print "total number of pseudogenes = $num_tot_pseudo_genes\n";
    print "\tunclassified pseudogenes = $num_unclass_pseudo_genes\n";
    print "\tprocessed pseudogenes = $num_proc_pseudo_genes\n";
    print "\tunprocessed pseudogenes = $num_unproc_pseudo_genes\n";
    print "Ig Pseudogenes = $num_pred_Ig_pseudogenes\n";
    print "Ig Segments = $num_Ig_segments\n\n\n";

    my @attribs;

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Known Genes',
     -CODE => 'KnownGeneCount',
     -VALUE => $num_known_genes,
     -DESCRIPTION => 'Total Number of Known Genes');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Novel CDS',
     -CODE => 'NovelCDSCount',
     -VALUE => $num_nov_CDS,
     -DESCRIPTION => 'Total Number of Novel CDSs');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Novel transcripts',
     -CODE => 'NovelTransCount',
     -VALUE => $num_nov_trans,
     -DESCRIPTION => 'Total Number of Novel Transcripts');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Putative transcripts',
     -CODE => 'PutTransCount',
     -VALUE => $num_put_trans,
     -DESCRIPTION => 'Total Number of Putative Transcripts');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Predicted transcripts',
     -CODE => 'PredTransCount',
     -VALUE => $num_pred_trans,
     -DESCRIPTION => 'Total Number of Predicted Transcripts');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Total Pseudogenes',
     -CODE => 'TotPsCount',
     -VALUE => $num_tot_pseudo_genes,
     -DESCRIPTION => 'Total Number of Pseudogenes');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Unclassified Pseudogenes',
     -CODE => 'UnclassPsCount',
     -VALUE => $num_unclass_pseudo_genes,
     -DESCRIPTION => 'Number of Unclassified Pseudogenes');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Processed Pseudogenes',
     -CODE => 'ProcPsCount',
     -VALUE => $num_proc_pseudo_genes,
     -DESCRIPTION => 'Number of Processed Pseudogenes');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Unprocessed Pseudogenes',
     -CODE => 'UnprocPsCount',
     -VALUE => $num_unproc_pseudo_genes,
     -DESCRIPTION => 'Number of Unprocessed Pseudogenes');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Ig Segments',
     -CODE => 'IgSegCount',
     -VALUE => $num_Ig_segments,
     -DESCRIPTION => 'Total Number of Ig Segments');

    push @attribs, Bio::EnsEMBL::Attribute->new
    (-NAME => 'Ig Pseudogene Segments',
     -CODE => 'IgPsSegCount',
     -VALUE => $num_pred_Ig_pseudogenes,
     -DESCRIPTION => 'Total Number of Ig Pseudogene Segments');

    if( $snps_present ) {
	my $snps = $slice->get_all_SNPs();
	push @attribs, Bio::EnsEMBL::Attribute->new
	(-NAME => 'SNP Count',
	 -CODE => 'SNPCount',
	 -VALUE => scalar( @$snps ),
	 -DESCRIPTION => 'Total Number of SNPs');
    }
    
    $attrib_adaptor->store_on_Slice($slice, \@attribs);
    #  print_chromo_stats([$slice]);
}



sub print_chromo_stats {
    my $chromosomes = shift;
    
    foreach my $chr (@$chromosomes) {
	print "\nchromosome: ",$chr->seq_region_name(),"\n";
	foreach my $attrib (@{$chr->get_all_Attributes()}) {
	    print "  ", $attrib->name(), ": ", $attrib->value(), "\n";
	}
    }
}


#
# tries to attach lite.
#

sub lite_attach {
    my $db = shift;
    
    my $core_db_name;
    $core_db_name = $db->dbname();
    if( $core_db_name !~ /_core_/ ) {
	return 0;
    }
    #
    # get a lost of all databases on that server
    #
    my $sth = $db->prepare( "show databases" );
    $sth->execute();
    my $all_db_names = $sth->fetchall_arrayref();
    my %all_db_names = map {( $_->[0] , 1)} @$all_db_names;
    my $snp_db_name = $core_db_name;
    $snp_db_name =~ s/_core_/_lite_/;
    if( ! exists $all_db_names{ $snp_db_name } ) {
	return 0;
    }
    
    my $snp_db = Bio::EnsEMBL::Lite::DBAdaptor->new
    ( -host => $db->host(),
      -user => $db->username(),
      -pass => $db->password(),
      -port => $db->port(),
      -dbname => $snp_db_name );
    $db->add_db_adaptor( "lite", $snp_db );
    return 1;
}


1;



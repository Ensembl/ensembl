#!/usr/local/lib/perl


my $locator = "Bio::EnsEMBL::DBSQL::Obj/host=obi-wan;dbname=ensembl;user=ensro;";
$db = Bio::EnsEMBL::DBLoader->new($locator);
@clones = $db->get_all_Clone_id();


foreach my $clone_id ( @clones ) {

    eval {
	my $clone = $db->get_Clone($clone_id);
	@genes = $clone->get_all_Genes();
	$genenumber = scalar @genes;
	my @exons;

	foreach $g ( @genes ) {
	    push(@exons,$g->each_unique_Exon);
	}

	$exonnumber= scalar @exons;
	$length = 0;
	foreach $exon ( @exons ) {
	    $length += $exon->lenght;
	}

	$gene{$clone_id} = $genenumber;
	$exon{$clone_id} = $exonnumber;
	$exonl{$clone_id} = $length;
    };

    if( $i > 20 ) {
	last;
    }

    $i++;
}

@clones = sort { $exon{$a} <=> $exon{$b} } keys %exon;

foreach $clone_id (@clones) {
    print "$clone_id\t$exon{$clone_id}\t$gene{$clone_id}\t$exonl{$clone_id}\n";
}



    

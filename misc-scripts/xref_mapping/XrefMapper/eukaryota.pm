
package XrefMapper::eukaryota;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


sub get_set_lists {

  return [["ExonerateGappedBest1", ["eukaryota","*"]]];

}

sub transcript_display_xref_sources {
    my $self     = shift;
    my $fullmode = shift;

    print STDERR "getting the list of external_dbs for assigning gene names from eukaryota.pm\n";

    my @list = qw(
                 RFAM
                 RNAMMER
                 TRNASCAN_SE
                 Uniprot_genename
                 ENA_GENE
                 BROAD_U_maydis
                 BROAD_F_oxysporum
                 BROAD_G_zeae
                 BROAD_G_moniliformis
                 BROAD_P_infestans
                 phyra_jgi_v1.1
                 physo1_jgi_v1.1
	   	 phatr_jgi_v2
		 phatr_jgi_v2_bd
                 PGD_GENE
                 Mycgr3_jgi_v2.0_gene
                 BROAD_Magnaporthe_DB
               );
    
    my %ignore;

    # Both methods
    
    if(!$fullmode){
	$ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';
    }
    else{
	$ignore{"EntrezGene"} = 'select ox.object_xref_id from object_xref ox, dependent_xref dx, source s1, xref x1, source s2, xref x2 where ox.object_xref_id = dx.object_xref_id and dx.dependent_xref_id = x1.xref_id and x1.source_id = s1.source_id and s1.name = "EntrezGene" and x2.xref_id = dx.master_xref_id and x2.source_id = s2.source_id and (s2.name like "Refseq_dna_predicted" or s2.name like "RefSeq_peptide_predicted") and ox.ox_status = "DUMP_OUT"';
	
    }
    
    return [\@list,\%ignore];
}

sub gene_description_sources {
  return (
          "PomBase_GENE",
          "PomBase_TRANSCRIPT",
          "Uniprot/SWISSPROT",
          "Uniprot/SPTREMBL",
          "BROAD_U_maydis",
          "BROAD_F_oxysporum",
          "BROAD_G_zeae",
          "BROAD_G_moniliformis",
          "BROAD_P_infestans",
          "phyra_jgi_v1.1",
          "physo1_jgi_v1.1",
	  "phatr_jgi_v2",
	  "phatr_jgi_v2_bd",
          "PGD_GENE",
          "BROAD_Magnaporthe_DB",
          "RFAM",
          "TRNASCAN_SE",
          "RNAMMER",
         );
}


sub set_source_id_to_external_name {
    
    my $self = shift;
    my $name_to_external_name_href = shift;

    my $source_id_to_external_name_href = {};
    my $name_to_source_id_href = {};
    
    print STDERR "overwritting set_source_id_to_external_name in eukaryota.pm\n";
    
    my $sql = 'select s.source_id, s.name from source s, xref x where x.source_id = s.source_id group by s.source_id'; # only get those of interest
    
    my $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    my ($id, $name);
    $sth->bind_columns(\$id, \$name);
    while($sth->fetch()){
	if(defined($name_to_external_name_href->{$name})){
	    $source_id_to_external_name_href->{$id} = $name;
	    $name_to_source_id_href->{$name} = $id;
	}
	elsif($name =~ /notransfer$/){
	}
	else{
	    die "ERROR: Could not find $name in external_db table please add this too continue";
	}
    }
    
    $sth->finish;
    
    return ($source_id_to_external_name_href, $name_to_source_id_href);
}

1;

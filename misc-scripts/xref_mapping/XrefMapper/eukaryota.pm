
package XrefMapper::eukaryota;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1';
  my %override_method_for_source = (  
	   ExonerateGappedBest5 => ['RefSeq_mRNA','RefSeq_mRNA_predicted', 'RefSeq_ncRNA', 'RefSeq_ncRNA_predicted' ],
         );

  return $default_method, \%override_method_for_source;
}

sub gene_display_xref_sources {
    my $self     = shift;

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
                 PGSC_GENE
                 PHYTOZOME_GMAX_GENE
               );
    
    my %ignore;

    
    #don't use EntrezGene labels dependent on predicted RefSeqs

    $ignore{'EntrezGene'} =<<IEG;
SELECT DISTINCT ox.object_xref_id
  FROM object_xref ox, dependent_xref dx, 
       xref xmas, xref xdep, 
       source smas, source sdep
    WHERE ox.xref_id = dx.dependent_xref_id AND
          dx.dependent_xref_id = xdep.xref_id AND
          dx.master_xref_id = xmas.xref_id AND
          xmas.source_id = smas.source_id AND
          xdep.source_id = sdep.source_id AND
          smas.name like "Refseq%predicted" AND
          sdep.name like "EntrezGene" AND
          ox.ox_status = "DUMP_OUT" 	 
IEG

    #don't use labels starting with LOC

    $ignore{'LOC_prefix'} =<<LOCP;
SELECT object_xref_id
  FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
   WHERE ox_status = 'DUMP_OUT' AND label REGEXP '^LOC[[:digit:]]+'
LOCP

    return [\@list,\%ignore];
}


sub transcript_display_xref_sources {
    my $self     = shift;

    print STDERR "getting the list of external_dbs for assigning transcript names from eukaryota.pm\n";

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
                 PGSC_GENE
                 PHYTOZOME_GMAX_GENE
               );
    
    my %ignore;

    
    #don't use EntrezGene labels dependent on predicted RefSeqs

    $ignore{'EntrezGene'} =<<IEG;
SELECT DISTINCT ox.object_xref_id
  FROM object_xref ox, dependent_xref dx, 
       xref xmas, xref xdep, 
       source smas, source sdep
    WHERE ox.xref_id = dx.dependent_xref_id AND
          dx.dependent_xref_id = xdep.xref_id AND
          dx.master_xref_id = xmas.xref_id AND
          xmas.source_id = smas.source_id AND
          xdep.source_id = sdep.source_id AND
          smas.name like "Refseq%predicted" AND
          sdep.name like "EntrezGene" AND
          ox.ox_status = "DUMP_OUT" 	 
IEG

    #don't use labels starting with LOC

    $ignore{'LOC_prefix'} =<<LOCP;
SELECT object_xref_id
  FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
   WHERE ox_status = 'DUMP_OUT' AND label REGEXP '^LOC[[:digit:]]+'
LOCP

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
          "PGSC_GENE",
          "PHYTOZOME_GMAX_GENE",
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

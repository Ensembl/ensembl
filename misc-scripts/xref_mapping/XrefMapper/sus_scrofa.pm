package XrefMapper::sus_scrofa;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists{

  return [["ExonerateGappedBest5", ["sus_scrofa","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["sus_scrofa","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["sus_scrofa","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["sus_scrofa","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["sus_scrofa","*"]]];

}


sub gene_description_sources {

  return ("RFAM",
          "RNAMMER",
          "TRNASCAN_SE",
	  "miRBase",
          "HGNC",
          "IMGT/GENE_DB",
	  "Uniprot/SWISSPROT",
	  "RefSeq_peptide",
	  "RefSeq_dna",
	  "Uniprot/Varsplic",
	  "Uniprot/SPTREMBL");

}

sub gene_description_filter_regexps {

  return ('^UNKNOWN\s+.*',
          '^undefined.*');

}
        

sub transcript_display_xref_sources {
  my $self     = shift;

  my @list = qw(HGNC
                MGI
                ZFIN_ID
                Clone_based_vega_gene
                Clone_based_ensembl_gene
                HGNC_transcript_name
                MGI_transcript_name
                ZFIN_ID_transcript_name
                Clone_based_vega_transcript
                Clone_based_ensembl_transcript
                miRBase
                RFAM
                RNAMMER
                IMGT/GENE_DB
                SGD
                flybase_symbol
                Anopheles_symbol
                Genoscope_annotated_gene
                Uniprot_genename
                Uniprot/SWISSPROT
                EntrezGene
                TRNASCAN_SE
                RefSeq_peptide
                RefSeq_mRNA
                RefSeq_ncRNA);

  my %ignore;

  $ignore{"EntrezGene"} =(<<'IEG');
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

  $ignore{"Uniprot/SPTREMBL"} =(<<BIGN);
SELECT object_xref_id
    FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
     WHERE ox_status = 'DUMP_OUT' AND name = 'Uniprot/SPTREMBL' 
      AND priority_description = 'protein_evidence_gt_2'
BIGN

  $ignore{"LOC"} =(<<CIGN);
SELECT object_xref_id
    FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
     WHERE ox_status = 'DUMP_OUT' AND name in ('Uniprot_genename','EntrezGene') 
      AND label like 'LOC%'
CIGN

  $ignore{"SSC."} =(<<DIGN);
SELECT object_xref_id
    FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
     WHERE ox_status = 'DUMP_OUT' AND name = 'Uniprot_genename' 
      AND label like 'SSC\.%'
DIGN
  return [\@list,\%ignore];
}

1;

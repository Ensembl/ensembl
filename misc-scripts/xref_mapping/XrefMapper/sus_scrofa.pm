package XrefMapper::sus_scrofa;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists{

  return [["ExonerateGappedBest5", ["xenopus_tropicalis","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["xenopus_tropicalis","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["xenopus_tropicalis","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["xenopus_tropicalis","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["xenopus_tropicalis","*"]]];

}

sub gene_description_filter_regexps {

  return ('^UNKNOWN\s+.*',
          '^undefined.*');

}



sub transcript_display_xref_sources {
  my $self     = shift;

  my @list = qw(HGNC
                MGI
                Clone_based_vega_gene
                Clone_based_ensembl_gene
                HGNC_transcript_name
                MGI_transcript_name
                Clone_based_vega_transcript
                Clone_based_ensembl_transcript
                miRBase
                RFAM
                IMGT/GENE_DB
                SGD
                flybase_symbol
                Anopheles_symbol
                Genoscope_annotated_gene
                Uniprot_genename
                Uniprot/SWISSPROT
                Uniprot/Varsplic
                Uniprot/SPTREMBL
                EntrezGene
                RefSeq_mRNA
                RefSeq_ncRNA
                RefSeq_dna
                RefSeq_peptide);

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

  $ignore{"EntrezGene/LOC"} =(<<BIGN);
SELECT object_xref_id
    FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
     WHERE ox_status = 'DUMP_OUT' AND name = 'EntrezGene' 
      AND accession = '%LOC%'
BIGN


  return [\@list,\%ignore];
}

1;

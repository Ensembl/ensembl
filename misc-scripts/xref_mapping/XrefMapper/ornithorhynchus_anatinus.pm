package XrefMapper::ornithorhynchus_anatinus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

# Same as in BasicMapper but Genoscope order reversed.

sub transcript_display_xref_sources {
  my @list = qw(Platypus_olfactory_receptor
		Oxford_FGU_Oa_tscript
		Oxford_FGU_Oa_gene
		RFAM
		miRBase
		IMGT/GENE_DB
		HUGO
		SGD
		MGI
		flybase_symbol
		Anopheles_symbol
		Genoscope_annotated_gene
		Uniprot/SWISSPROT
		Uniprot/Varsplic
		RefSeq_peptide
		RefSeq_dna
		Uniprot/SPTREMBL
		EntrezGene);

  my %ignore;
  $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';
  
  return [\@list,\%ignore];
}

sub gene_description_sources {

 return ("RFAM",
         "miRBase",
         "IMGT/GENE_DB",
         "Uniprot/SWISSPROT",
         "RefSeq_peptide",
         "RefSeq_dna",
         "Uniprot/Varsplic",
         "Uniprot/SPTREMBL");

}




1;

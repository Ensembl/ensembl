package XrefMapper::ornithorhynchus_anatinus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

# Same as in BasicMapper but Genoscope order reversed.

sub transcript_display_xref_sources {
  print "RETURNING THE display xref sources\n";
  return ('Platypus_olfactory_receptor',
          'Oxford_FGU_Oa_tscript',
          'Oxford_FGU_Oa_gene',
          'RFAM',
          'miRBase',
          'IMGT/GENE_DB',
          'HUGO',
          'SGD',
          'MarkerSymbol',
          'flybase_symbol',
          'Anopheles_symbol',
          'Genoscope_annotated_gene',
          'Uniprot/SWISSPROT',
          'Uniprot/Varsplic',
          'RefSeq_peptide',
          'RefSeq_dna',
          'Uniprot/SPTREMBL',
          'EntrezGene');

}

1;

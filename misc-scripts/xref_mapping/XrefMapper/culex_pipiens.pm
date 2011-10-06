
package XrefMapper::culex_pipiens;

use  XrefMapper::BasicMapper;
use  XrefMapper::VBCoordinateMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1_culex", ["culex_pipiens","*"]]];

}

# transcript, gene display_xrefs can use defaults
# since anopheles_symbol is "before" Uniprot

# If there is an Anopheles_symbol xref, use its description
sub gene_description_sources {

  return ("Anopheles_symbol",
	  "Uniprot/SWISSPROT",
	  "RefSeq_peptide",
	  "RefSeq_dna",
	  "Uniprot/SPTREMBL",
	  #"RefSeq_peptide_predicted",
	  #"RefSeq_dna_predicted",
	  #"EntrezGene");
          );
}

sub transcript_display_xref_sources {

  my @list = qw(RFAM
		miRBase
		Uniprot/SWISSPROT
		Uniprot/Varsplic
		Uniprot/SPTREMBL);

  my %ignore;

    $ignore{"Uniprot/SPTREMBL"} =(<<BIGN);
SELECT object_xref_id
    FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
     WHERE ox_status = 'DUMP_OUT' AND name = 'Uniprot/SPTREMBL' 
      AND priority_description = 'protein_evidence_gt_3'
BIGN

  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {

  return ();

}

1;

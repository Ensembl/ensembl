package XrefMapper::tetraodon_nigroviridis;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

# Same as in BasicMapper but Genoscope order reversed.

sub transcript_display_xref_sources {

  my @list = qw(HGNC
		MGI
		wormbase_transcript
		flybase_symbol
		Anopheles_symbol
		Genoscope_annotated_gene
		Genoscope_predicted_transcript
		Uniprot/SWISSPROT
		RefSeq
		Uniprot/SPTREMBL
		LocusLink);

  my %ignore;

  $ignore{"Uniprot/SPTREMBL"} =(<<BIGN);
SELECT object_xref_id
    FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
     WHERE ox_status = 'DUMP_OUT' AND name = 'Uniprot/SPTREMBL' 
      AND priority_description = 'protein_evidence_gt_2'
BIGN


  return [\@list,\%ignore];


}

1;

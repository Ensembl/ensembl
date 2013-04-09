package XrefMapper::tetraodon_nigroviridis;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub gene_display_xref_sources {
  my $self     = shift;
	
  my @list = qw(Genoscope_annotated_gene
                RFAM
                miRBase
                Uniprot_genename
                EntrezGene
                Genoscope_pred_gene
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

  my @list = qw(	
	        Genoscope_ann_transcript
		RFAM
                miRBase          
		Uniprot/SWISSPROT
                Uniprot/Varsplic
                Genoscope_pred_transcript
);

  my %ignore;

  return [\@list,\%ignore];

}

1;

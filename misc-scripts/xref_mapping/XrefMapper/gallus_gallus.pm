package XrefMapper::gallus_gallus;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };
use strict;


sub gene_display_xref_sources {
  my $self     = shift;

  my @list = qw(CGNC
                RFAM
                miRBase
                Uniprot_gn
                EntrezGene);

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
          ox.ox_status = "DUMP_OUT" AND
          ox.master_xref_id = dx.master_xref_id 
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

  return ("miRBase",
	  "RFAM", 
	  "CGNC",
	  "Uniprot/SWISSPROT", 
	  "Uniprot/Varsplic", 
	  "RefSeq_peptide", 
	  "RefSeq_mRNA");

}


1;

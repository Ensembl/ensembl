package XrefMapper::aedes_aegypti;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  #return [ ["ExonerateGappedBest1", ["aedes_aegypti","*"]] , ["ExonerateGappedBest1_aedesGnBk", ["aedes_aegypti","AedesGenBank"]] ];  #Initial + new AedesGenBank - doesn't work
  #return [ ["ExonerateGappedBest1", ["aedes_aegypti","Uniprot/SWISSPROT"]] ,
  #	   ["ExonerateGappedBest1", ["aedes_aegypti","Uniprot/SPTREMBL"]] ,
  #	   ["ExonerateGappedBest1", ["aedes_aegypti","UniGene"]] ,
  #	   ["ExonerateGappedBest1", ["aedes_aegypti","EMBL"]] ,
  #	   ["ExonerateGappedBest1", ["aedes_aegypti","PDB"]] ,
  #	   ["ExonerateGappedBest1", ["aedes_aegypti","protein_id"]] ,
  #	   ["ExonerateGappedBest1", ["aedes_aegypti","GO"]] ,
  #	   ["ExonerateGappedBest1", ["aedes_aegypti","Interpro"]] ,
  #	   ["ExonerateGappedBest1_aedesGnBk", ["aedes_aegypti","AedesGenBank"]] ];  #Initial + new AedesGenBank - doesn't work
  return [["ExonerateGappedBest1", ["aedes_aegypti","*"]]];                        #Initial - OK
  #return [["ExonerateGappedBest1_aedesGnBk", ["aedes_aegypti","AedesGenBank"]]];    #New AedesGnBk - OK
}


sub gene_description_sources {
  return ("Uniprot/SWISSPROT",
	  "RefSeq_peptide",
	  "RefSeq_dna",
	  "Uniprot/SPTREMBL");
}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {
  return ();
}1;

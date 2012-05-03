package XrefMapper::aedes_aegypti;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };



sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1';
  my %override_method_for_source = (  
	   ExonerateGappedBest5 => ['RefSeq_mRNA','RefSeq_mRNA_predicted', 'RefSeq_ncRNA', 'RefSeq_ncRNA_predicted' ],
         );

  return $default_method, \%override_method_for_source;
}



#Reverse order: last one has higher priority!
sub gene_description_sources {
  return ("VB_External_Description",
	  "Uniprot/SWISSPROT",
          "VB_Community_Annotation"
	  );
}

sub transcript_display_xref_sources {
  my @list = qw(
		Uniprot/SWISSPROT
		VB_Community_Annotation
                );

  my %ignore;
  #$ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';

  return [\@list,\%ignore];

}


1;

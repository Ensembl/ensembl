package XrefMapper::anopheles_gambiae;

use  XrefMapper::BasicMapper;
use  XrefMapper::VBCoordinateMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };



sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1_55_perc_id';
  my %override_method_for_source = (  
	   ExonerateGappedBest5_55_perc_id => ['RefSeq_mRNA','RefSeq_mRNA_predicted', 'RefSeq_ncRNA', 'RefSeq_ncRNA_predicted' ],
         );

  return $default_method, \%override_method_for_source;
}



# transcript, gene display_xrefs can use defaults
# since anopheles_symbol is "before" Uniprot

# If there is an Anopheles_symbol xref, use its description

# mh4 says Anopheles_symbol doesn't get chosen over UniP
# (but they do get chosen in other cases)

sub gene_description_sources {

  return ("VB_Community_Annotation",
	  "Uniprot/SWISSPROT",
	  "VB_RNA_Description",
          );
}

sub gene_display_xref_sources {

  my @list = qw(RFAM
		miRBase
		VB_Community_Annotation
		Uniprot/SWISSPROT
		VB_RNA_Description
	       );

  my %ignore;
  return [\@list,\%ignore];

}

sub transcript_display_xref_sources {

  my @list = qw(RFAM
		miRBase
		VB_Community_Annotation
		Uniprot/SWISSPROT
		VB_RNA_Description
	       );

  my %ignore;
  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {

  return ();

}

1;

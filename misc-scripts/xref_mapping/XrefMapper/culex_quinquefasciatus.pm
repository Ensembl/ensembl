package XrefMapper::culex_quinquefasciatus;

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

#Reverse order: the latest one is the list has higher precedence!

sub gene_description_sources {

  return (
	  "VB_External_Description",
	  "VB_RNA_Description",
	  "Uniprot/SWISSPROT",
	  "VB_Community_Annotation"
          );
}

sub transcript_display_xref_sources {

  my @list = qw(RFAM
		miRBase
		Uniprot/SWISSPROT
		VB_Community_Annotation
             );

  my %ignore;
  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {

  return ();

}

sub no_source_label_list{
  my $self = shift;
  my @list;

  print "Using no_source_label_list :-)\n";
  #foreach my $ex (qw("VB RNA Description" "VB External Description")){
  #  $list{$ex} = 1;
  #}	

  push @list,"VB_RNA_Description";
  push @list,"VB_External_Description";

  return \@list;
}

1;

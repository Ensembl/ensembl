=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

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



#Reverse order: last one has higher priority! <- This comment is wrong (mikkel 12/10/2012)
# index [0] has higher priority! 
sub gene_description_sources {
  return qw(
  	     VB_Community_Annotation
  	     Uniprot/SWISSPROT
  	     VB_External_Description
	     );
}
#Old order VB_External_Description Uniprot/SWISSPROT VB_Community_Annotation
# only sources in this list are used when setting  display_xref_id
sub transcript_display_xref_sources {
  my @list = qw(
  	  	VB_Community_Annotation
  	  	Uniprot/SWISSPROT
  	  	VB_External_Description
		);
#old VB_External_Description Uniprot/SWISSPROT VB_Community_Annotation
  my %ignore;
  #$ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';

  return [\@list,\%ignore];

}


1;

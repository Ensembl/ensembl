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

package XrefMapper::ixodes_scapularis;

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

#Reverse order: second one has higher priority!
sub gene_description_sources {
  return qw(
  	     VB_Community_Annotation
  	     Uniprot/SWISSPROT
  	     VB_External_Description
	     );
}

#Reverse order: second one has higher priority!
sub transcript_display_xref_sources {
  my $self     = shift;
  my $fullmode = shift;

  my @list = qw(
  	  	VB_Community_Annotation
  	  	Uniprot/SWISSPROT
  	  	VB_External_Description
		);

  my %ignore;

  return [\@list,\%ignore];

}

# regexps to match any descriptons we want to filter out
sub gene_description_filter_regexps {
  return ();
}1;

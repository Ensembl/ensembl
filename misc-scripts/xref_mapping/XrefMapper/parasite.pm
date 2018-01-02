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

package XrefMapper::parasite;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


# This module is activated by specifying "taxon=wormbase" in the mapping input file
# It contains some common config for worms maintained by WormBase (i.e. having genes
# with WBGene ids etc)

sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1';
  my %override_method_for_source = (  
           ExonerateGappedBest5 => ['RefSeq_mRNA',
                                    'RefSeq_mRNA_predicted', 
                                    'RefSeq_ncRNA', 
                                    'RefSeq_ncRNA_predicted' ],
         );

  return $default_method, \%override_method_for_source;
}


sub transcript_names_from_gene {
  return;
}


sub gene_description_sources {

  return ("RFAM",
          "RNAMMER",
          "TRNASCAN_SE",
	  "miRBase",
          "HGNC",
          "IMGT/GENE_DB",
	  "Uniprot/SWISSPROT",
	  "RefSeq_peptide",
          "Uniprot/SPTREMBL",
      );

}


sub gene_description_filter_regexps {

  return (
           '^Uncharacterized protein\s*',
           '^Putative uncharacterized protein\s*',
           '^Hypothetical protein\s*',
   );

}

1;

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

package XrefMapper::saccharomyces_cerevisiae;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1';
  my %override_method_for_source = (  
	   ExonerateGappedBest5 => ['RefSeq_mRNA','RefSeq_mRNA_predicted', 'RefSeq_ncRNA', 'RefSeq_ncRNA_predicted' ],
         );

  return $default_method, \%override_method_for_source;
}


# Cerevisiae is imported from SGD. The gene and transcript stable IDs
# are the SGD identifiers. The display_xref_ids for genes and
# transcripts are calculated directly rather than via the more complex
# priority-based method in BasicMapper.pm

sub build_display_xrefs {

  my ($self, $type, $external_db) = @_;

  print "Setting $type display_xrefs from $type stable IDs\n";
  my $dir = $self->core()->dir();

  my $sql = "UPDATE $type t, xref x, external_db e SET t.display_xref_id=x.xref_id WHERE t.stable_id=x.dbprimary_acc AND e.external_db_id=x.external_db_id AND e.db_name=\'${external_db}\'\n";

  open (SQL, ">$dir/${type}_display_xref.sql");

  print SQL $sql;

  close(SQL);

}

sub gene_display_xref_sources {
    my $self     = shift;

    my @list = qw(
                SGD_GENE
               );
    
    my %ignore;
     
    return [\@list,\%ignore];
}


sub transcript_display_xref_sources {
    my $self     = shift;

    my @list = qw(
                SGD_TRANSCRIPT
               );
    
    my %ignore;
     
    return [\@list,\%ignore];
}

sub gene_description_sources {
  return (
          "SGD_GENE"
         );
}


sub gene_description_filter_regexps {

  return ();

}



1;

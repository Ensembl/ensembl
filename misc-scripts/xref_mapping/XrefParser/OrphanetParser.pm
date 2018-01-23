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

package XrefParser::OrphanetParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use XML::LibXML;

use base qw( XrefParser::BaseParser );

sub run {


  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my $xml_file = @{$files}[0];
  
  print STDERR "Orphanet file to parse, $xml_file\n" if($verbose);

  my %gene_disorders; 

  my $term = undef;
  my $desc = undef;
  
  my $xml_parser = XML::LibXML->new();
  my $orphanet_doc    = $xml_parser->parse_file($xml_file);

  my ($jdbor_node) = $orphanet_doc->findnodes('JDBOR');
  my $release = $jdbor_node->getAttribute('version');
  # Set release
  $self->set_release( $source_id,$release );


  my $gene_disorder_count = 0;
  my $disorder_count = 0;
  
  foreach my $disorder ($orphanet_doc->findnodes('JDBOR/DisorderList/Disorder')) {
    my ($orpha_number_node) = $disorder->findnodes('./OrphaNumber');
    my $orpha_number = $orpha_number_node->to_literal;
    my ($name_node) = $disorder->findnodes('./Name');
    my $name = $name_node->to_literal;

    my @genes = $disorder->findnodes('./DisorderGeneAssociationList/DisorderGeneAssociation/Gene');

    if ( scalar(@genes) > 0) {

	$disorder_count++;
    }

    foreach my $gene (@genes) {

	my $ref;
	#get the HGNC xref
	foreach my $external_reference_node ($gene->findnodes('./ExternalReferenceList/ExternalReference')) {
	    my ($source_node) = $external_reference_node->findnodes('./Source');
	    if ($source_node->to_literal =~ /HGNC/) {
		my ($ref_node) = $external_reference_node->findnodes('./Reference');
		$ref = $ref_node->to_literal;
	    }
	}
	if (defined($ref)) {
	    $gene_disorders{$ref}{$orpha_number} = $name;
	    $gene_disorder_count++;
	}
    }
    
  }
 
  print  "Parsed $disorder_count disorders\n";
  print  "Found $gene_disorder_count genes associated with disorders\n";


  #get the mapping that are already there so that we don't get lots of duplicates.
  # stored in the global hash xref_dependent_mapped.
  $self->get_dependent_mappings($source_id);


  my (%hgnc) = %{$self->get_valid_codes("HGNC", $species_id)};

  print "got " . scalar(keys (%hgnc)) . " HGNC entries\n";
   
  print "species_id, source_id: $species_id, $source_id\n";

  my $added = 0;
  my %Orphanet_xrefs;
 
  foreach my $hgnc_acc (keys (%gene_disorders)) {
    
      # Get the master_xref_id
    
      if(!defined($hgnc{$hgnc_acc})){
	print STDERR "failed to get the master_xref_if for HGNC, $hgnc_acc!\n";
      }
      else{
	foreach my $master_xref_id (@{$hgnc{$hgnc_acc}}){
	  
	  foreach my $orpha_number (keys %{$gene_disorders{$hgnc_acc}}) {

	    my $description = $gene_disorders{$hgnc_acc}{$orpha_number};    

	    $self->add_dependent_xref({ master_xref_id => $master_xref_id,
					acc            => $orpha_number,
					label          => $description || $orpha_number,
					source_id      => $source_id,
					species_id     => $species_id,
				        desc           => $description });
	    $Orphanet_xrefs{$orpha_number}++;
	    $added++;
	  }
	}
      }
    
  }

  print "Added " . scalar(keys %Orphanet_xrefs) . " Orphanet xrefs and " . $added . " dependent xrefs\n";
  if ($added > 0) {
      return 0;
  } else {
      return 1;
  }

}

1;

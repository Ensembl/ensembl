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

package XrefParser::PHIbaseParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
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

  my $phi_xml_file = @{$files}[0];
  
  print STDERR "PhiBase file to parse, $phi_xml_file\n" if($verbose);

  my %phi_mapping;
  my %taxIds;

  my $term = undef;
  my $desc = undef;
  
  my $phi_parser = XML::LibXML->new();
  my $phi_doc    = $phi_parser->parse_file($phi_xml_file);

  my $concepts_count = 0;
  
  foreach my $concept ($phi_doc->findnodes('ondex:ondexdata/ondexdataseq/concepts/concept')) {
    my ($pid_node) = $concept->findnodes('./pid');
    my $pid = $pid_node->to_literal;
    my $concept_accessions_aref = $concept->findnodes('./coaccessions/concept_accession');
    my $uniprot_acc = undef;

    foreach my $concept_accession (@$concept_accessions_aref) {
	
	# get the one associated with UniProt
	my ($elementOf) = $concept_accession->findnodes('./elementOf');

	if ($elementOf->to_literal =~ /UPROT/) {
	    my ($accession) = $concept_accession->findnodes('./accession');
	    $uniprot_acc = $accession->to_literal;
	}
    }

    if (!defined $uniprot_acc) {
	# print STDERR "phi id, $pid, no uniprot mapping found in xml file!\n";
    }
    else {
      if (!defined $phi_mapping{$uniprot_acc}) {
	$phi_mapping{$uniprot_acc} = 
	  {
	   -phi_ids => [$pid],
	   -tax_id  => undef,
	  };
      }
      else {
	my $phi_href = $phi_mapping{$uniprot_acc};
	my $aref = $phi_href->{-phi_ids};
	push (@$aref, $pid);
      }
    }

      # Get the TaxId

    my $concept_gds_aref = $concept->findnodes('./cogds/concept_gds');
    my $taxId = undef;

    foreach my $concept_gds (@$concept_gds_aref) {

      # get the one associated with Taxid
      my ($attrname) = $concept_gds->findnodes('./attrname');
      if ($attrname->to_literal =~ /TAXID/) {
	my ($value) = $concept_gds->findnodes('./value');
	$taxId = $value->to_literal;
	$taxId =~ s/\D//g;
	
	if (! defined $taxIds{$taxId}) {
	  $taxIds{$taxId} = 1;
	}
	
	if (defined $uniprot_acc) {
	  my $phi_href = $phi_mapping{$uniprot_acc};
	  $phi_href->{-tax_id} = $taxId;
	}
      }
    }
    $concepts_count++;
  }

  my @phis = keys (%phi_mapping);

  print  "Parsed $concepts_count concepts\n";
  print  "Found " . @phis . " with UniProt mapping!\n";

  print  "Found " . keys (%taxIds) . " different taxIds\n";

  #get the "main" PHIbase source id.
  $source_id = $self->get_source_id_for_source_name("PHIbase");


  #get the mapping that are already there so that we don't get lots of duplicates.
  # stored in the global hash xref_dependent_mapped.
  $self->get_dependent_mappings($source_id);

  #if(!defined($species_id)){
  #  $species_id = $self->get_species_id_for_filename($phi_xml_file);
  #}

  my $swiss_miss=0;
  my (%swiss) = %{$self->get_valid_codes("uniprot/", $species_id)};

  print "got " . keys (%swiss) . " Uniprot entries\n";
   
  print "species_id, source_id: $species_id, $source_id\n";

  # Don't check only the species_id, but all taxIds specified in xref_config.ini

  my %species2tax = $self->species_id2taxonomy();
  my @tax_ids = @{$species2tax{$species_id}};

  print "tax_ids from xref_config.ini file: " . join (', ', @tax_ids) . "\n";

  my $added = 0;

  foreach my $uniprot_acc (keys (%phi_mapping)) {
    my $phis_href = $phi_mapping{$uniprot_acc};
    my $taxId = $phis_href->{-tax_id};
    if (grep {$_ eq $taxId} @tax_ids) {
      # Get the master_xref_id
      # and the linkage

      if(!defined($swiss{$uniprot_acc})){
	print STDERR "failed to get the master_xref_if for UniProt, $uniprot_acc!\n";
	# one reason it happens is that the UniProt identifier is attached to a different tax node in Phibase that it is in UniProt
      }
      else{
	foreach my $master_xref_id (@{$swiss{$uniprot_acc}}){
	  #	      print STDERR "master_xref_id, $master_xref_id\n";
	  my $linkage = undef;
	  my $phis_aref = $phis_href->{-phi_ids};
	  foreach my $phibase_id (@$phis_aref) {
	    #		  print STDERR "Adding xrefs for phibase id, $phibase_id\n";
	    $self->add_dependent_xref({ master_xref_id => $master_xref_id,
					acc            => $phibase_id,
					label          => $phibase_id,
					source_id      => $source_id,
					species_id     => $species_id });
	    $added++;
	  }
	}
      }
    }
  }

  print "Added $added PHIbase xrefs\n";

  return 0;

}

1;

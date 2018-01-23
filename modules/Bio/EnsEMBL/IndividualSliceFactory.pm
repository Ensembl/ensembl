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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::IndividualSliceFactory;

use strict;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Scalar::Util qw(weaken);

=head2 new
=cut

sub new{
    my $caller = shift;
    my $class = ref($caller) || $caller;

    #creates many IndividualSlice objects from the Population

    my ($population_name, $coord_system, $start, $end, $strand, $seq_region_name, $seq_region_length, $adaptor) = rearrange(['POPULATION', 'COORD_SYSTEM','START','END','STRAND','SEQ_REGION_NAME','SEQ_REGION_LENGTH', 'ADAPTOR'],@_);

    my $self = bless {
	population_name => $population_name,
	coord_system => $coord_system,
	start => $start,
	end => $end,
	strand => $strand,
	seq_region_name => $seq_region_name,
	seq_region_length => $seq_region_length},$class;

    $self->adaptor($adaptor);
    return $self;
}

sub adaptor {
  my $self = shift;

  if(@_) {
    my $ad = shift;
    if($ad && (!ref($ad) || !$ad->isa('Bio::EnsEMBL::DBSQL::BaseAdaptor'))) {
      throw('Adaptor argument must be a Bio::EnsEMBL::DBSQL::BaseAdaptor');
    }
    weaken($self->{'adaptor'} = $ad);
  }

  return $self->{'adaptor'}
}

sub get_all_IndividualSlice{
    my $self = shift;

    my $slice;
    if(!$self->adaptor) {
	warning('Cannot get IndividualSlice features without attached adaptor');
	return '';
    }
    my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

    unless($variation_db) {
	warning("Variation database must be attached to core database to " .
		"retrieve variation information" );
	return '';
    }
    #get the AlleleFeatures in the Population
    my $af_adaptor = $variation_db->get_AlleleFeatureAdaptor;
    
    if( $af_adaptor ) {
	#set the adaptor to retrieve data from genotype table
	$af_adaptor->from_IndividualSlice(1);
	#get the Individual for the given strain
	my $population_adaptor = $variation_db->get_PopulationAdaptor;
	my $individual_adaptor = $variation_db->get_IndividualAdaptor;
	if ($population_adaptor && $individual_adaptor){
	    $slice = Bio::EnsEMBL::Slice->new(-coord_system => $self->{'coord_system'},
					      -start => $self->{'start'},
					      -end => $self->{'end'},
					      -strand => $self->{'strand'},
					      -seq_region_name => $self->{'seq_region_name'},
					      -seq_region_length => $self->{'seq_region_length'},
					      -adaptor => $self->adaptor
					      );
	    my $population = $population_adaptor->fetch_by_name($self->{'population_name'}); 
	    #check that there is such population in the database
	    if (defined $population){
		#get all the AlleleFeatures in the $population and the Slice given
		my $allele_features = $af_adaptor->fetch_all_by_Slice($slice,$population);
		#get Individuals in the Population
		my $individuals = $individual_adaptor->fetch_all_by_Population($population);		
		return $self->_rearrange_Individuals_Alleles($individuals,$allele_features);
	    }
	    else{ 
		warning("Population not in the database");
		return '';
 
	    }
	}
	else{
	    warning("Not possible to retrieve PopulationAdaptor from the variation database");
	    return '';
	}	
    }
    
    else{
	warning("Not possible to retrieve AlleleFeatureAdaptor from variation database");
	return '';
    }
}

sub _rearrange_Individuals_Alleles{
    my $self = shift;
    my $individuals = shift;
    my $allele_features;
    my $individual_slice;
    #create the hash with all the individuals
    my %individuals_ids;
    #foreach of the individual, create the IndividualSlice object and add it to the mapping hash
    foreach my $individual (@{$individuals}){
	$individual_slice = Bio::EnsEMBL::Variation::IndividualSlice->new(
	    -coord_system => $self->{'coord_system'},
	    -start => $self->{'$start'},
	    -end  => $self->{'end'},
	    -strand => $self->{'strand'},
	    -seq_region_name => $self->{'seq_region_name'},
	    -seq_region_length => $self->{'seq_region_length'},
	    -individual => $individual->name);
	
	$individuals_ids{$individual->dbID} = $individual_slice;
    }

    #and rearrange all the AlleleFeatures to the individuals
    foreach my $allele_feature (@{$allele_features}){
	$individuals_ids{$allele_feature->{'_sample_id'}}->add_AlleleFeature($allele_feature);
    }
    my @result = values %individuals_ids;
    return \@result;
}


1;

package Bio::EnsEMBL::IndividualSliceFactory;

use strict;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::IndividualSlice;
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Time::HiRes qw(gettimeofday tv_interval);
use Data::Dumper;

=head2 new
=cut

sub new{
    my $caller = shift;
    my $class = ref($caller) || $caller;

    #creates many IndividualSlice objects from the Population

    my ($coord_system, $start, $end, $strand, $seq_region_name, $seq_region_length, $adaptor) = rearrange(['COORD_SYSTEM','START','END','STRAND','SEQ_REGION_NAME','SEQ_REGION_LENGTH', 'ADAPTOR'],@_);

    return bless {
	coord_system => $coord_system,
	start => $start,
	end => $end,
	strand => $strand,
	seq_region_name => $seq_region_name,
	seq_region_length => $seq_region_length,
	adaptor => $adaptor},$class;
}

sub get_all_IndividualSlice{
    my $self = shift;

    my $slice;
    if(!$self->{'adaptor'}) {
	warning('Cannot get IndividualSlice features without attached adaptor');
	return '';
    }
    my $variation_db = $self->{'adaptor'}->db->get_db_adaptor('variation');

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
	my $individual_adaptor = $variation_db->get_IndividualAdaptor;
	if ($individual_adaptor){
	    $slice = Bio::EnsEMBL::Slice->new(-coord_system => $self->{'coord_system'},
					      -start => $self->{'start'},
					      -end => $self->{'end'},
					      -strand => $self->{'strand'},
					      -seq_region_name => $self->{'seq_region_name'},
					      -seq_region_length => $self->{'seq_region_length'},
					      -adaptor => $self->{'adaptor'}
					      );
	    #get all the AlleleFeatures in the $population and the Slice given
	    my $allele_features = $af_adaptor->fetch_all_IndividualSlice($slice);

	    return $self->_rearrange_Individuals_Alleles($allele_features);	    
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
    my $allele_features = shift;
    my $individual_slice;
    #create the hash with all the individuals
    my %individuals_ids;
    my $individual;
    my $variation_db = $self->{'adaptor'}->db->get_db_adaptor('variation');
    my $individual_adaptor = $variation_db->get_IndividualAdaptor();
    #and rearrange all the AlleleFeatures to the individuals
    foreach my $allele_feature (@{$allele_features}){
	if (!defined $individuals_ids{$allele_feature->{'_sample_id'}}){
	    $individual = $individual_adaptor->fetch_by_dbID($allele_feature->{'_sample_id'});

	    #create the IndividualSlice
	    $individual_slice = Bio::EnsEMBL::IndividualSlice->new(
								   -coord_system => $self->{'coord_system'},
								   -start => $self->{'start'},
								   -end  => $self->{'end'},
								   -strand => $self->{'strand'},
								   -seq_region_name => $self->{'seq_region_name'},
								   -seq_region_length => $self->{'seq_region_length'},
								   -individual => $individual->name,
								   -sample_id => $individual->dbID);
	    $individuals_ids{$allele_feature->{'_sample_id'}} = $individual_slice;
	}
	$individuals_ids{$allele_feature->{'_sample_id'}}->add_AlleleFeature($allele_feature);
    }
    my @result = values %individuals_ids;
    return \@result;
}


1;

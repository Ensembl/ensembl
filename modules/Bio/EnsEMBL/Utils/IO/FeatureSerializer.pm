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

=pod


=head1 NAME

FeatureSerializer - An abstract serializer for turning EnsEMBL Features into other formats

=head1 SYNOPSIS

my $serializer = new Bio::EnsEMBL::Utils::IO::FeatureSerializer( $filehandle );
$serializer->print_feature_list( \@list_of_features );

=head1 DESCRIPTION

Generic class for serializing features Subclass this class to create a 
format-specific serializer. Be sure to implement print_feature at the 
bare minimum

=cut

package Bio::EnsEMBL::Utils::IO::FeatureSerializer;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Utils::IO::Serializer);

# Majority of methods inherited from Serializer

=head2 print_feature

	Arg []     : None
	Description: Abstract method to print a feature.
                     Implemented in derived classes.
        Returntype : None

=cut

sub print_feature {
	throw( "print_feature method not implemented.");
}

=head2 print_feature_list

	Arg [1]    : Listref of Features
	Description: Run print_feature on every feature in the list
        Returntype : None
=cut

sub print_feature_list {
	my $self = shift;
	my $feature_list = shift;
	if (ref($feature_list) eq 'ARRAY') {
		if (scalar(@$feature_list) > 0) {$self->{'achieved_something'} = 1;} 
		foreach my $feature (@{$feature_list}) {
			$self->print_feature($feature);
		}
	}
	else {
		throw( "print_feature_list requires a listref as argument" );
	}
}

=head2 print_feature_Iterator

	Arg [1]    : Bio::EnsEMBL::Utils::Iterator
	Description: Automatically spools through an iterator for convenience
	Returntype : None
=cut

sub print_feature_Iterator {
	my $self = shift;
	my $iterator = shift;
	if ($iterator->can('has_next')) {
		$iterator->each(
			sub {
				$self->print_feature($_); 
				$self->{'achieved_something'} = 1;
			}
		);
	}
	else {
		throw("Supplied iterator does not look like Bio::EnsEMBL::Utils::Iterator");
	}
}

1;

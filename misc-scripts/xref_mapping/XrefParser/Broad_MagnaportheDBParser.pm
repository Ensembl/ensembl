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

package XrefParser::Broad_MagnaportheDBParser;

use strict;
use warnings;
use Carp;

use base qw( XrefParser::Broad_GenericParser );

sub run {

    my ($self, $ref_arg) = @_;
    my $source_id    = $ref_arg->{source_id};
    my $species_id   = $ref_arg->{species_id};
    my $files        = $ref_arg->{files};
    my $verbose      = $ref_arg->{verbose};
    
    if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
	croak "Need to pass source_id, species_id and  files as pairs";
    }
    $verbose |= 0;
    
    my $file = @{$files}[0];
    
    my $broad_source_id = $self->get_source_id_for_source_name("BROAD_Magnaporthe_DB");
    $ref_arg->{source_id} = $broad_source_id;

    # Parse the Broad gene summary file

    my $xrefs_aref = $self->parse ($ref_arg);

    # Store the parsed data as direct xrefs

    my $errorcode = $self->store ($xrefs_aref);

    return $errorcode;

}

1;


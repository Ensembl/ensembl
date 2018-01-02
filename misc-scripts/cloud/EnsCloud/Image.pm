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

package EnsCloud::Image;
use Moose;


has volume_bundle => (isa => 'EnsCloud::Image::VolumeBundle', is=>'rw', required => 1);

has 'root_device_snapshot' => (isa => 'Str', is => 'rw') ;
has 'root_device' => (isa => 'Str', is => 'rw') ;

has 'kernel' => (isa => 'Str', is => 'rw') ;
has 'species' => (isa => 'Str', is => 'rw', required => 1) ;
has 'tag' => (isa => 'Str', is => 'rw', required => 1) ;

1;

__END__

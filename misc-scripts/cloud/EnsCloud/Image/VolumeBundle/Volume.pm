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

package EnsCloud::Image::VolumeBundle::Volume;
use Moose;
use Data::Dumper;
has 'tag' => ( isa => 'Str', is => 'rw' );
has 'species' => ( isa => 'Str', is => 'rw', required => 1 );


has 'total_size' => (
    traits => ['Number'],
    is     => 'rw',

    #      isa     => 'Int',
    default => 0,
    handles => {
        add_size => 'add',
    }
);

has 'volume_id'   => ( isa => 'Str', is => 'rw' );
has 'ensembl_release'   => ( isa => 'Int', is => 'rw' );
has 'snapshot_id' => ( isa => 'Str', is => 'rw' );
has 'status' => ( isa => 'Str', is => 'rw', required =>1 , default => 'no snapshot' );
has 'device'      => ( isa => 'Str', is => 'rw' );
has 'dbs'         => (
    traits  => ['Array'],
    handles => {
        all_dbs => 'elements',
        add_db  => 'push',
        next_db => 'shift',
	sort_in_place_curried => [sort_in_place => sub {$_[0]->type cmp $_[1]->type}]
    },
    isa => 'ArrayRef[EnsCloud::Image::VolumeBundle::Volume::DatabaseDetails]',
);

1;

__END__

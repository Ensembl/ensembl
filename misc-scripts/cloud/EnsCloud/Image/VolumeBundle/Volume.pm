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

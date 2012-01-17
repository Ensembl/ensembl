package EnsCloud::Image::VolumeBundle;

use Moose;

has 'volumes' => (
    traits  => ['Array'],
    handles => {
        all_volumes => 'elements',
        add_volumes => 'push',
        next_volume => 'shift',
    },
    isa => 'ArrayRef[EnsCloud::Image::VolumeBundle::Volume]',

);

1;

__END__

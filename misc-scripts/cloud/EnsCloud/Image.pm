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

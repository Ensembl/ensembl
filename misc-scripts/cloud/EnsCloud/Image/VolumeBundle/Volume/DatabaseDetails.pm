package EnsCloud::Image::VolumeBundle::Volume::DatabaseDetails;
use Moose;
use CHI;


has 'myd_path' => (isa => 'Str' , is => 'rw', required => 1);
has 'name' => (isa => 'Str', is => 'rw', required => 1 ) ;
has 'size' => (isa => 'Str', is => 'rw', required => 1);
has 'type' => (isa => 'Str', is => 'rw', required => 1);
has 'copy_path' => (isa => 'Str', is => 'rw', );
has 'is_copied' => (isa => 'Str', is => 'rw', required => 1, default => 0);








1;

__END__

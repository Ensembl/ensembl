package Bio::EnsEMBL::DBSQL::MapSet;
use strict;
use vars qw( $AUTOLOAD );

## Instantiates object of class Bio::EnsEMBL::DBSQL::MapSet;
sub new {
    my $class = shift;
    my $self = {
        'mapset_id'     => $_[0],   'name'          => $_[1],
        'code'          => $_[2],   'description'   => $_[3]
    };
    bless $self,$class;
    return $self;
}
 
## Let the ID be the name of the mapset
sub id {
    my $self = shift;
    return $self->{'name'};
}

sub destroy { return 1; }

sub AUTOLOAD {
    my $self = shift;
    no strict 'refs';
    my $var = $AUTOLOAD;
    $var =~ s/.*:://;
    return $self->{$var}
}

1;
package Bio::EnsEMBL::DBSQL::MapFrag;
use strict;
use vars qw( $AUTOLOAD );

sub new {
    my $class = shift;
    my $self = {
        '_synonyms'   => [],        '_annotations'=> {},
        '_mapsets'    => {},
        '_vc_start'   => $_[0],     'mapfrag_id'  => $_[1],     
        'type'        => $_[2],     'seq'         => $_[3],
        'seq_type'    => $_[4],     'seq_start'   => $_[5],
        'seq_end'     => $_[6],     'orientation' => $_[7],
        'name'        => $_[8]
    };
    bless $self,$class;
    return $self;
}
 
sub add_synonym {
    my( $self, $value ) = @_;
    push @{$self->{'_synonyms'}}, $value;
}

sub add_mapset {
    my( $self, $value ) = @_;
    $self->{'_mapsets'}{$value->{'code'}} = $value;
}

sub add_annotation {
    my( $self, $type, $value ) = @_;
    $self->{'_annotations'}{$type} = $value;
}

sub synonyms {
    my $self = shift;
    return @{$self->{'_synonyms'}};
}

sub mapsets {
    my $self = shift;
    return @{$self->{'_mapsets'}};
}

sub is_in_mapset {
    my( $self, $mapset ) = @_;
    return exists( $self->{'_mapsets'}{$mapset} ) ? 1 : 0;
}

sub length {
    my $self = shift;
    return $self->{'seq_end'} - $self->{'seq_start'} + 1;
}

sub start { # these are in VC co-ordinates
    my $self = shift;
    return $self->{'seq_start'} - $self->{'_vc_start'} + 1;
}

sub end {   # these are in VC co-ordinates
    my $self = shift;
    return $self->{'seq_end'} - $self->{'_vc_start'} + 1;
}

sub id {
    my $self = shift;
    return $self->name;
}

sub destroy { return 1; }

sub bacinfo {
    my $self = shift;
    return ('',
            'Interpolated',
            'Both ends located',
            'Start located',
            'End located')
        [$self->{'_annotations'}{'BACend_flag'}];
}

sub AUTOLOAD {
    my $self = shift;
    no strict 'refs';
    my $var = $AUTOLOAD;
    $var =~ s/.*:://;
    return exists $self->{$var} ? $self->{$var} : $self->{'_annotations'}{$var};
}

1;
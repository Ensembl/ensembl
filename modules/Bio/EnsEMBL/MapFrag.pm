package Bio::EnsEMBL::MapFrag;
use strict;
use vars qw( $AUTOLOAD );

## Instantiates object of class Bio::EnsEMBL::MapFrag;

sub new {
    my $class = shift;
    my $self = {
        '_synonyms'   => {},        '_embl_accs'  => {},
        '_annotations'=> {},        '_mapsets'    => {},        
        '_vc_start'   => $_[0],     'mapfrag_id'  => $_[1],     
        'type'        => $_[2],     'seq'         => $_[3],
        'seq_type'    => $_[4],     'seq_start'   => $_[5],
        'seq_end'     => $_[6],     'orientation' => $_[7],
        'name'        => $_[8]
    };
    bless $self,$class;
    return $self;
}

## Functions for adding entries to the synonyms, embl_accs,
##           annotations and mapsets hashses
sub add_synonym {
    my( $self, $value ) = @_;
    $self->{'_synonyms'}{$value} = 1;
}

sub add_embl_acc {
    my( $self, $value ) = @_;
    $self->{'_embl_accs'}{$value} = 1;
}

sub add_mapset {
    my( $self, $mapset ) = @_;
    $self->{'_mapsets'}{ $mapset->code } = $mapset;
}

sub add_annotation {
    my( $self, $type, $value ) = @_;
    $self->{'_annotations'}{$type} = $value;
}

## Retrieve information from the synonmys, embl_accs and mapsets hashes.
sub synonyms {
    my $self = shift;
    return keys %{$self->{'_synonyms'}};
}

sub embl_accs {
    my $self = shift;
    return keys %{$self->{'_embl_accs'}};
}
sub mapsets {
    my $self = shift;
    return values %{$self->{'_mapsets'}};
}

## Check to see if the clone is in the given mapset
sub is_in_mapset {
    my( $self, $mapset ) = @_;
    return exists( $self->{'_mapsets'}{$mapset} ) ? 1 : 0;
}

## Return the lenghth of the "clone"
sub length {
    my $self = shift;
    return $self->{'seq_end'}   - $self->{'seq_start'} + 1;
}

## Return the start/end of the clone in VC co-oridnates
sub start {
    my $self = shift;
    return $self->{'seq_start'} - $self->{'_vc_start'} + 1;
}

sub end {   # these are in VC co-ordinates
    my $self = shift;
    return $self->{'seq_end'}   - $self->{'_vc_start'} + 1;
}

## Let the ID be the name of the clone
sub id {
    my $self = shift;
    return $self->{'name'};
}

sub destroy { return 1; }

## Returns a textual description of what the bacflags variable means
sub bacinfo {
    my $self = shift;
    return ('Interpolated',
            'Start located',
            'End located',
            'Both ends located')
        [$self->{'_annotations'}{'BACend_flag'}];
}

## Deals with autgenerating accessors for other variables and elements of
##           annotations hash.
sub AUTOLOAD {
    my $self = shift;
    no strict 'refs';
    (my $var = $AUTOLOAD) =~ s/.*:://;
    if(exists $self->{$var}) {
        *{$AUTOLOAD} = sub { return $_[0]{$var}; };
        return $self->{$var}
    } elsif($self->{'_annotations'}{$var}) {
        *{$AUTOLOAD} = sub { return $_[0]{'_annotations'}{$var}; };
        return $self->{'_annotations'}{$var}
    }
}

1;

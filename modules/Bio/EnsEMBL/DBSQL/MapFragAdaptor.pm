package Bio::EnsEMBL::DBSQL::MapFragAdaptor;
use strict;
use vars '@ISA';
use Bio::EnsEMBL::DBSQL::MapFrag;
use Bio::EnsEMBL::DBSQL::MapSet;

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

## Overide new class so that we can create the 
## cache and max_feature_length values.........
sub new {
    my $class = shift;
    my $self = $class->SUPER::new( @_ );
    $self->max_feature_length( 1e8 );
    $self->{'_cache'} = {};
    return $self;
}

sub _get_dnafrag_id { 
    my($self,$chr_name) = @_;
    
    my $sth = $self->prepare("select dnafrag_id from dnafrag where name = ?");
    $sth->execute($chr_name);
    my ($dnafrag_id) = $sth->fetchrow_array();

    return $dnafrag_id;
}

sub _get_mapset_id {
    my($self,$mapset_code) = @_;

    my $sth = $self->prepare("select mapset_id from mapset where code = ?");
    $sth->execute($mapset_code);
    my ($mapset_id) = $sth->fetchrow_array();

    return $mapset_id;
}

sub fetch_mapset_chr_start_end {
    my( $self, $mapset_code, $chr_name, $chr_start, $chr_end ) = @_;
    my $key = join ':', $mapset_code, $chr_name, $chr_start, $chr_end;
    
    return @{$self->{'_cache'}{$key}} if $self->{'_cache'}{$key};

    my( $dnafrag_id, $mapset_id );
    
    unless(  ( $dnafrag_id = $self->_get_dnafrag_id( $chr_name ) ) &&
             ( $mapset_id  = $self->_get_mapset_id( $mapset_code ) )  ) {
        $self->{'_cache'}{$key} = [];
        return ();
    }

    my $sth = $self->prepare(
        qq( select mf.mapfrag_id, mf.type, mf.name,
                   mf.seq_start, mf.seq_end, mf.orientation,
                   df.name as seq, df.dnafrag_type as seq_type,
                   mat.code as note_type, ma.value as note
              from mapfrag as mf, dnafrag as df, mapfrag_mapset as mm
                   left join mapannotation as ma on ma.mapfrag_id = mf.mapfrag_id
                   left join mapannotationtype as mat on ma.mapannotationtype_id = mat.mapannotationtype_id
             where mm.mapset_id = ? and mm.mapfrag_id = mf.mapfrag_id and mf.dnafrag_id = df.dnafrag_id
                    ).
            ( $chr_name ? "and df.dnafrag_id = ?".( $chr_start ? " and mf.seq_start <= ? and mf.seq_start >= ? and mf.seq_end >= ?" : "" ) : "").
        qq(  order by mf.mapfrag_id, mat.code )
    );
        
    $sth->execute(
        $mapset_id, (
             $dnafrag_id ?
            ( $dnafrag_id, ( $chr_start ?
                           ($chr_end, $chr_start - $self->max_feature_length, $chr_start) :
                           () ) ) :
            ()
        )
    );
    my @map_frags = ();
    my $old_id    = 0;
    my $map_frag  = undef;
    while( my $data = $sth->fetchrow_hashref() ) {
    if($data->{'mapfrag_id'}!=$old_id) {
        push @map_frags, $map_frag if defined $map_frag;
            $map_frag = Bio::EnsEMBL::DBSQL::MapFrag->new(
                $chr_start || 1,
                $data->{'mapfrag_id'},      
                $data->{'type'},            $data->{'seq'},
                $data->{'seq_type'},        $data->{'seq_start'},
                $data->{'seq_end'},         $data->{'orientation'},
                $data->{'name'},
            );
            $old_id = $data->{'mapfrag_id'}
        }
        if($data->{'note_type'} eq 'synonym') {
            $map_frag->add_synonym( $data->{'note'} ); 
        } elsif($data->{'note_type'} eq 'embl_acc') {
            $map_frag->add_embl_acc( $data->{'note'} ); 
        } else {
            $map_frag->add_annotation( $data->{'note_type'}, $data->{'note'} );
        }
    }
    push @map_frags, $map_frag if defined $map_frag;
    $self->{'_cache'}{$key} = \@map_frags;
    return @map_frags;
}

sub fetch_by_internal_id {
    my( $self, $ID) = @_;
    my $key = "ID:$ID";
    print STDERR "FETCH_BY INTERNAL ID $ID\n";
    return $self->{'_cache'}{$key} if $self->{'_cache'}{$key};
    my $sth = $self->prepare(
        qq( select mf.mapfrag_id, mf.type, mf.name,
                   mf.seq_start, mf.seq_end, mf.orientation,
                   df.name as seq, df.dnafrag_type as seq_type,
                   mat.code as note_type, ma.value as note
              from mapfrag as mf, dnafrag as df
                   left join mapannotation as ma on ma.mapfrag_id = mf.mapfrag_id
                   left join mapannotationtype as mat on ma.mapannotationtype_id = mat.mapannotationtype_id
             where mf.dnafrag_id = df.dnafrag_id and mf.mapfrag_id = ? )
    );
        
    $sth->execute( $ID );
    return undef unless($sth->rows>0);
    
    my $map_frag  = undef;
    while( my $data = $sth->fetchrow_hashref() ) {
        unless($map_frag) {
            $map_frag = Bio::EnsEMBL::DBSQL::MapFrag->new(
                1,
                $data->{'mapfrag_id'},      
                $data->{'type'},            $data->{'seq'},
                $data->{'seq_type'},        $data->{'seq_start'},
                $data->{'seq_end'},         $data->{'orientation'},
                $data->{'name'},
            );
        }
        if($data->{'note_type'} eq 'synonym') {
            $map_frag->add_synonym( $data->{'note'} ); 
        } elsif($data->{'note_type'} eq 'embl_acc') {
            $map_frag->add_embl_acc( $data->{'note'} ); 
        } else {
            $map_frag->add_annotation( $data->{'note_type'}, $data->{'note'} );
        }
    }

    $sth = $self->prepare(
        qq(select ms.id, ms.code, ms.name, ms.description
             from mapset as ms, mapfrag_mapset as mm
            where ms.mapset_id = mm.mapset_id and mm.mapfrag_id = ?
        )
    );
    $sth->execute( $ID );

    while( my $data = $sth->fetchrow_arrayref() ) {
        $map_frag->add_mapset( Bio::EnsEMBL::DBSQL::MapSet->new( $data ) )
    }
    
    return $self->{'_cache'}{$key} = $map_frag;
}

sub fetch_by_embl_acc {
    my( $self, $embl_acc ) = @_;
    my $key = "name:$embl_acc";
    return $self->{'_cache'}{$key} if $self->{'_cache'}{$key};

    my $sth = $self->prepare(
        "select ma.mapfrag_id
           from mapannotation as ma,
                mapannotationtype as mat
          where ma.value = ? and mat.code = 'embl_acc' and
                ma.mapannotationtype_id = mat.mapannotationtype_id"
    );
    $sth->execute( $embl_acc );
    my( $ID ) = $sth->fetchrow_array();
    return $self->{'_cache'}{$key} = $self->fetch_by_internal_id( $ID );    
}

sub fetch_by_name {
    my( $self, $name) = @_;
    print STDERR "FETCHING BY NAME\n";
    my $key = "name:$name";
    return $self->{'_cache'}{$key} if $self->{'_cache'}{$key};

    my $sth = $self->prepare(
        "select mapfrag_id
           from mapfrag
          where name = ?"
    );
    $sth->execute( $name );
    my( $ID ) = $sth->fetchrow_array();
    return $self->{'_cache'}{$key} = $self->fetch_by_internal_id( $ID );    
}

sub fetch_by_synonym {
    my( $self, $synonym ) = @_;
    print STDERR "FETCHING BY SYNONYM\n";
    my $key = "name:$synonym";
    return $self->{'_cache'}{$key} if $self->{'_cache'}{$key};
    my $sth = $self->prepare(
        "select ma.mapfrag_id
           from mapannotation as ma,
                mapannotationtype as mat
          where ma.value = ? and mat.code in ('embl_acc', 'synonym') and
                ma.mapannotationtype_id = mat.mapannotationtype_id"
    );
    $sth->execute( $synonym );
    my( $ID ) = $sth->fetchrow_array();
    if($ID) {
        return $self->{'_cache'}{$key} = $self->fetch_by_internal_id( $ID );    
    } else {
        return $self->fetch_by_name( $synonym );    
    }
}

sub max_feature_length {
    my $self = shift;
    $self->{'_max_feature_length'} = shift if( @_ );
    return $self->{'_max_feature_length'} || 1e9;
}

sub get_mapsets {
    my $self = shift;
    my $flag = shift;
    my $sth = $self->prepare(
        "select mapset_id, code, name, description
           from mapset"
    );
    $sth->execute();
    my %results = ();
    my $key = $flag eq 'mapset_id' ? 'mapset_id' : 'code'; # Key to store hash on...
    while( my $data = $sth->fetchrow_arrayref() ) {
        $results{ $data->{$key} } = Bio::EnsEMBL::DBSQL::MapSet->new( $data );
    }
    return %results;
}

sub has_mapset {
    my $self = shift;
    my $name = shift;
    my $sth = $self->prepare(
        "select 1 from mapset as ms, mapfrag_mapset as mm where ms.code = ? and ms.mapset_id = mm.mapset_id limit 1"
    );
    $sth->execute( $name );
    return $sth->fetchrow_array();
}
1;
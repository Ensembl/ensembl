package Bio::EnsEMBL::DBSQL::MapFragAdaptor;
use strict;
use vars '@ISA';
use Bio::EnsEMBL::MapFrag;
use Bio::EnsEMBL::MapSet;

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

## Overide new class so that we can create the 
## cache and max_feature_length values.........


=head2 new

  Arg [1]    : list @args 
               superclass constructor arguments
  Example    : none
  Description: Creates a new MapFragAdaptor.  Superclass constructor is
               overridden so that max_feature_length attribute can be 
               initialized.
  Returntype : Bio::EnsEMBL::DBSQL::MapFragAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBConnection

=cut

sub new {
    my $class = shift;
    my $self = $class->SUPER::new( @_ );
    $self->max_feature_length( 1e9 );
    $self->{'_cache'} = {};
    return $self;
}


=head2 _get_dnafrag_id

  Arg [1]    : string $chr_name
  Example    : none
  Description: PRIVATE method
               returns a dnafrag id via its name 
  Returntype : int
  Exceptions : none
  Caller     : internal

=cut

sub _get_dnafrag_id { 
    my($self,$chr_name) = @_;
    
    my $sth = $self->prepare("select dnafrag_id from dnafrag where name = ?");
    $sth->execute($chr_name);
    my ($dnafrag_id) = $sth->fetchrow_array();

    return $dnafrag_id;
}


=head2 _get_mapset_id

  Arg [1]    : string $mapset_code
  Example    : none
  Description: PRIVATE retrieves a mapset id via its code
  Returntype : int
  Exceptions : none
  Caller     : internal

=cut

sub _get_mapset_id {
    my($self,$mapset_code) = @_;

    my $sth = $self->prepare("select mapset_id from mapset where code = ?");
    $sth->execute($mapset_code);
    my ($mapset_id) = $sth->fetchrow_array();

    return $mapset_id;
}



=head2 fetch_by_mapset_chr_start_end

  Arg [1]    : string $mapset_code
               The code of the mapset to retrieve map frags from
  Arg [2]    : string $chr_name (optional)
               The name of the chromosome or dna_fragment to obtain map
               frags from
  Arg [3]    : int $chr_start (optional)
               The start of the region to obtain map_frags from
  Arg [4]    : int $chr_end (optional)
               The end of the region to obtain map_frags from
  Example    : @mfs = $mf_adaptor->fetch_by_mapset_chr_start_end('Tilepath');
  Description: Retrieves a list of MapFragments from a given mapset within an
               optionally specified region of the assembly 
  Returntype : Bio::EnsEMBL::MapFrag
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice

=cut

sub fetch_by_mapset_chr_start_end {
    my( $self, $mapset_code, $chr_name, $chr_start, $chr_end ) = @_;
    my $key = join ':', $mapset_code, $chr_name, $chr_start, $chr_end;
    
    return @{$self->{'_cache'}{$key}} if $self->{'_cache'}{$key};

    my( $dnafrag_id, $mapset_id );
    
    unless( 
        ( $dnafrag_id = $self->_get_dnafrag_id( $chr_name ) ) &&
        ( $mapset_id  = $self->_get_mapset_id( $mapset_code ) )
    ) {
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
            ( $chr_name ? "and df.dnafrag_id = ?".( defined($chr_start) ? " and mf.seq_start <= ? and mf.seq_start >= ? and mf.seq_end >= ?" : "" ) : "").
        qq(  order by mf.mapfrag_id, mat.code )
    );
        
    $sth->execute(
        $mapset_id, (
             $dnafrag_id ?
            ( $dnafrag_id, ( defined($chr_start) ?
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
            $map_frag = Bio::EnsEMBL::MapFrag->new(
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


=head2 fetch_by_dbID

  Arg [1]    : int $ID
               the internal database id of the map frag to obtain 
  Example    : $map_frag = $map_frag_adaptor->fetch_by_dbID(123);
  Description: Obtains a map fragment object via its unique database identifier
  Returntype : Bio::EnsEMBL::MapFrag
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_dbID {
    my( $self, $ID) = @_;
    my $key = "ID:$ID";

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
            $map_frag = Bio::EnsEMBL::MapFrag->new(
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
        qq(select ms.mapset_id, ms.code, ms.name, ms.description
             from mapset as ms, mapfrag_mapset as mm
            where ms.mapset_id = mm.mapset_id and mm.mapfrag_id = ?
        )
    );
    $sth->execute( $ID );

    while( my $data = $sth->fetchrow_arrayref() ) {
        $map_frag->add_mapset( Bio::EnsEMBL::MapSet->new( $data ) )
    }
    
    return $self->{'_cache'}{$key} = $map_frag;
}


=head2 fetch_by_embl_acc

  Arg [1]    : string $embl_acc
               the embl accesion number
  Example    : $map_frag = $map_frag_adaptor->fetch_by_embl_acc(AC092813);
  Description: Retrieves a mapfragment via its associated EMBL accession
  Returntype : Bio::EnsEMBL::MapFrag
  Exceptions : none
  Caller     : general

=cut

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
    return $self->{'_cache'}{$key} = $self->fetch_by_dbID( $ID );    
}


=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $map_frag = $map_frag_adaptor->fetch_by_name('NT_032954');
  Description: Retrieves a map fragment via its name
  Returntype : Bio::EnsEMBL::MapFrag
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_name {
    my( $self, $name) = @_;

    my $key = "name:$name";
    return $self->{'_cache'}{$key} if $self->{'_cache'}{$key};

    my $sth = $self->prepare(
        "select mapfrag_id
           from mapfrag
          where name = ?"
    );
    $sth->execute( $name );
    my( $ID ) = $sth->fetchrow_array();
    return $self->{'_cache'}{$key} = $self->fetch_by_dbID( $ID );    
}


=head2 fetch_by_synonym

  Arg [1]    : string $synonym
               the synonym for the desirec map fragment
  Example    : $map_frag = $map_frag_adaptor->fetch_by_synonym('bA269M20'); 
  Description: Retrieves a map fragment via its synonym
  Returntype : Bio::EnsEMBL::MapFrag
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_synonym {
    my( $self, $synonym ) = @_;

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
        return $self->{'_cache'}{$key} = $self->fetch_by_dbID( $ID );    
    } else {
        return $self->fetch_by_name( $synonym );    
    }
}


=head2 max_feature_length

  Arg [1]    : int $length (optional)
               The new maximum feature length value
  Example    : $max_feature_length = $map_frag_adaptor->max_feature_length();
  Description: Getter/Setter for the maximum feature length which is used
               to limit the region of internal SQL queries for performance 
               reasons.  The default value is 1e9
  Returntype : int
  Exceptions : none
  Caller     : internal

=cut

sub max_feature_length {
    my $self = shift;
    $self->{'_max_feature_length'} = shift if( @_ );
    return $self->{'_max_feature_length'} || 1e9;
}


=head2 get_mapsets

  Arg [1]    : string $flag (optional) 
               if set to 'mapset_id' then the returned hash of data is 
               hashed on the mapset_id rather than on the code 
  Example    : %mapsets = $map_frag_adaptor->get_mapsets();
  Description: Retrieves a list of mapsets hashed on either their code 
               (default) or mapset_id. 
  Returntype : hash of Bio::EnsEMBL::MapSets
  Exceptions : none
  Caller     : general

=cut

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
        $results{ $data->{$key} } = Bio::EnsEMBL::MapSet->new( $data );
    }
    return %results;
}


=head2 has_mapset

  Arg [1]    : string $name
               the name of the mapset to check for
  Example    : if($map_frag_adaptor->has_mapset('cloneset')) do something; 
  Description: Returns true if a mapset with code $name exists
  Returntype : boolean
  Exceptions : none
  Caller     : general

=cut

sub has_mapset {
    my $self = shift;
    my $name = shift;
    my $sth = $self->prepare(
        "select 1 from mapset as ms, mapfrag_mapset as mm 
         where ms.code = ? and ms.mapset_id = mm.mapset_id limit 1"
    );
    $sth->execute( $name );
    return $sth->fetchrow_array();
}
1;

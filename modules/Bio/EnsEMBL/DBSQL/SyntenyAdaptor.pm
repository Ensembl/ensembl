package Bio::EnsEMBL::DBSQL::SyntenyAdaptor;
use vars qw( @ISA );
use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

## Sets the database name and the two species to compare

sub setSpecies {
    my( $self, $synteny_db, $species1, $species2 ) = @_;
    
    $self->{'_synteny_db'}        = $synteny_db;
    $self->{'_species_main'}      = $species1;
    $self->{'_species_secondary'} = $species2;
}

## Grabs the synteny blocks back.

sub get_synteny_for_chromosome {
    my( $self, $chr ) = @_; # if chr = undef return all synteny pairs
    
    my @data = ();
    my $SYNTENY_DB = $self->{'_synteny_db'};
    my $sth =$self->prepare(
        "select sr.synteny_region_id,
                df.dnafrag_type as core_type,  df.name as core_name, 
                dfr.seq_start as core_start,   dfr.seq_end as core_end, 
                df_h.dnafrag_type as hit_type, df_h.name as hit_name, 
                dfr_h.seq_start as hit_start,  dfr_h.seq_end as hit_end,
                sr.rel_orientation
           from ${SYNTENY_DB}.dnafrag as df,         ${SYNTENY_DB}.dnafrag as df_h, 
                ${SYNTENY_DB}.dnafrag_region as dfr, ${SYNTENY_DB}.dnafrag_region as dfr_h, 
                ${SYNTENY_DB}.genome_db as gd,       ${SYNTENY_DB}.genome_db as gd_h,
                ${SYNTENY_DB}.synteny_region as sr
          where gd.name = ?   and gd.genome_db_id   = df.genome_db_id and
                gd_h.name = ? and gd_h.genome_db_id = df_h.genome_db_id and
                df.dnafrag_id   = dfr.dnafrag_id and
                df_h.dnafrag_id = dfr_h.dnafrag_id and
                dfr.synteny_region_id   = sr.synteny_region_id and
                dfr_h.synteny_region_id = sr.synteny_region_id".
          ( defined $chr ? " and df.name = ? " : '' ).
          "order by df.name, dfr.seq_start          "
    );
    $sth->execute($self->{'_species_main'}, $self->{'_species_secondary'}, "$chr" );
    while(my $Q = $sth->fetchrow_arrayref()) {
        push @data, {
            'synteny_id'    => $Q->[0],
            'seq_type'      => $Q->[1],     'chr_name'      => $Q->[2],
            'chr_start'     => $Q->[3],     'chr_end'       => $Q->[4],
            'hit_seq_type'  => $Q->[5],     'hit_chr_name'  => $Q->[6],
            'hit_chr_start' => $Q->[7],     'hit_chr_end'   => $Q->[8],
            'rel_ori'       => $Q->[9]
        }
    }
    return \@data;
}    
